"""Property-based tests for benchmark.metrics.compute_error_metrics.

**Validates: Requirements 4.2, 4.3, 4.4, 4.5**

Property 3: Error metrics correctness
For any two equal-length value arrays and tolerance t > 0, verify max/mean
relative error and flagging logic.
"""

import math
from typing import Dict, List, Tuple

import pytest
from hypothesis import given, assume, settings
from hypothesis import strategies as st

from benchmark.metrics import compute_error_metrics
from benchmark.types import ErrorMetrics


# --- Strategies ---

# Strategy for a single finite, non-NaN float value (bounded to avoid overflow)
finite_float = st.floats(min_value=-1e100, max_value=1e100, allow_nan=False, allow_infinity=False)

# Strategy for a positive tolerance
positive_tolerance = st.floats(min_value=1e-15, max_value=1e10, allow_nan=False, allow_infinity=False)

# Strategy for a non-zero finite float
nonzero_float = st.floats(min_value=-1e100, max_value=1e100, allow_nan=False, allow_infinity=False).filter(
    lambda x: x != 0.0
)


def layer_value_dicts(
    min_layers: int = 1,
    max_layers: int = 5,
    min_species: int = 1,
    max_species: int = 5,
) -> st.SearchStrategy[Tuple[Dict[int, List[float]], Dict[int, List[float]]]]:
    """Generate a pair of layer->values dicts with matching keys and list lengths."""

    @st.composite
    def _build(draw):
        num_layers = draw(st.integers(min_value=min_layers, max_value=max_layers))
        num_species = draw(st.integers(min_value=min_species, max_value=max_species))
        layer_keys = draw(
            st.lists(
                st.integers(min_value=1, max_value=60),
                min_size=num_layers,
                max_size=num_layers,
                unique=True,
            )
        )

        fort_dict: Dict[int, List[float]] = {}
        cpp_dict: Dict[int, List[float]] = {}

        for key in layer_keys:
            fort_vals = draw(
                st.lists(finite_float, min_size=num_species, max_size=num_species)
            )
            cpp_vals = draw(
                st.lists(finite_float, min_size=num_species, max_size=num_species)
            )
            fort_dict[key] = fort_vals
            cpp_dict[key] = cpp_vals

        return fort_dict, cpp_dict

    return _build()


class TestErrorMetricsPropertyIdenticalInputs:
    """Property: For identical inputs, max and mean errors are 0 and flagged is False."""

    @given(
        data=st.data(),
        tolerance=positive_tolerance,
    )
    @settings(max_examples=100)
    def test_identical_inputs_zero_error(self, data, tolerance):
        """When fortran and cpp values are identical, errors must be 0.0 and flagged False."""
        num_layers = data.draw(st.integers(min_value=1, max_value=5))
        num_species = data.draw(st.integers(min_value=1, max_value=5))
        layer_keys = data.draw(
            st.lists(
                st.integers(min_value=1, max_value=60),
                min_size=num_layers,
                max_size=num_layers,
                unique=True,
            )
        )

        values: Dict[int, List[float]] = {}
        for key in layer_keys:
            vals = data.draw(
                st.lists(finite_float, min_size=num_species, max_size=num_species)
            )
            values[key] = vals

        result = compute_error_metrics(values, values, tolerance)

        assert result.max_relative_error == 0.0
        assert result.mean_relative_error == 0.0
        assert result.flagged is False


class TestErrorMetricsPropertyRelativeError:
    """Property: For non-zero reference, relative error = |cpp - fort| / |fort|."""

    @given(
        data=st.data(),
        tolerance=positive_tolerance,
    )
    @settings(max_examples=100)
    def test_relative_error_computation(self, data, tolerance):
        """Relative error for non-zero reference values is |cpp - fort| / |fort|."""
        num_species = data.draw(st.integers(min_value=1, max_value=4))
        # Generate a single layer with non-zero fortran values
        layer_key = data.draw(st.integers(min_value=1, max_value=60))

        fort_vals = data.draw(
            st.lists(nonzero_float, min_size=num_species, max_size=num_species)
        )
        cpp_vals = data.draw(
            st.lists(finite_float, min_size=num_species, max_size=num_species)
        )

        fort_dict = {layer_key: fort_vals}
        cpp_dict = {layer_key: cpp_vals}

        result = compute_error_metrics(fort_dict, cpp_dict, tolerance)

        # Compute expected errors manually
        expected_errors = [
            abs(cpp_vals[i] - fort_vals[i]) / abs(fort_vals[i])
            for i in range(num_species)
        ]
        expected_max = max(expected_errors)
        expected_mean = sum(expected_errors) / len(expected_errors)

        assert result.max_relative_error == pytest.approx(expected_max, rel=1e-10)
        assert result.mean_relative_error == pytest.approx(expected_mean, rel=1e-10)
        assert result.num_zero_reference == 0


class TestErrorMetricsPropertyZeroReference:
    """Property: For zero reference values, absolute difference is used."""

    @given(
        data=st.data(),
        tolerance=positive_tolerance,
    )
    @settings(max_examples=100)
    def test_zero_reference_uses_absolute_diff(self, data, tolerance):
        """When fortran value is 0, error = |cpp - fort| = |cpp|."""
        num_species = data.draw(st.integers(min_value=1, max_value=4))
        layer_key = data.draw(st.integers(min_value=1, max_value=60))

        # All fortran values are 0
        fort_vals = [0.0] * num_species
        cpp_vals = data.draw(
            st.lists(finite_float, min_size=num_species, max_size=num_species)
        )

        fort_dict = {layer_key: fort_vals}
        cpp_dict = {layer_key: cpp_vals}

        result = compute_error_metrics(fort_dict, cpp_dict, tolerance)

        # Expected: absolute difference |cpp - 0| = |cpp|
        expected_errors = [abs(cpp_vals[i]) for i in range(num_species)]
        expected_max = max(expected_errors)
        expected_mean = sum(expected_errors) / len(expected_errors)

        assert result.max_relative_error == pytest.approx(expected_max, rel=1e-10)
        assert result.mean_relative_error == pytest.approx(expected_mean, rel=1e-10)
        assert result.num_zero_reference == num_species


class TestErrorMetricsPropertyMaxError:
    """Property: max_relative_error is the maximum of all per-element errors."""

    @given(
        pairs=layer_value_dicts(),
        tolerance=positive_tolerance,
    )
    @settings(max_examples=100)
    def test_max_error_is_element_maximum(self, pairs, tolerance):
        """max_relative_error equals the max of individually computed per-element errors."""
        fort_dict, cpp_dict = pairs

        result = compute_error_metrics(fort_dict, cpp_dict, tolerance)

        # Manually compute all per-element errors
        errors: List[float] = []
        common_layers = sorted(set(fort_dict.keys()) & set(cpp_dict.keys()))
        for layer in common_layers:
            fort_vals = fort_dict[layer]
            cpp_vals = cpp_dict[layer]
            num_species = min(len(fort_vals), len(cpp_vals))
            for i in range(num_species):
                if fort_vals[i] == 0.0:
                    errors.append(abs(cpp_vals[i] - fort_vals[i]))
                else:
                    errors.append(abs(cpp_vals[i] - fort_vals[i]) / abs(fort_vals[i]))

        if errors:
            assert result.max_relative_error == pytest.approx(max(errors), rel=1e-10)
        else:
            assert result.max_relative_error == 0.0


class TestErrorMetricsPropertyMeanError:
    """Property: mean_relative_error is the arithmetic mean of all per-element errors."""

    @given(
        pairs=layer_value_dicts(),
        tolerance=positive_tolerance,
    )
    @settings(max_examples=100)
    def test_mean_error_is_arithmetic_mean(self, pairs, tolerance):
        """mean_relative_error equals the arithmetic mean of per-element errors."""
        fort_dict, cpp_dict = pairs

        result = compute_error_metrics(fort_dict, cpp_dict, tolerance)

        # Manually compute all per-element errors
        errors: List[float] = []
        common_layers = sorted(set(fort_dict.keys()) & set(cpp_dict.keys()))
        for layer in common_layers:
            fort_vals = fort_dict[layer]
            cpp_vals = cpp_dict[layer]
            num_species = min(len(fort_vals), len(cpp_vals))
            for i in range(num_species):
                if fort_vals[i] == 0.0:
                    errors.append(abs(cpp_vals[i] - fort_vals[i]))
                else:
                    errors.append(abs(cpp_vals[i] - fort_vals[i]) / abs(fort_vals[i]))

        if errors:
            expected_mean = sum(errors) / len(errors)
            assert result.mean_relative_error == pytest.approx(expected_mean, rel=1e-10)
        else:
            assert result.mean_relative_error == 0.0


class TestErrorMetricsPropertyFlagging:
    """Property: flagged is True iff max_relative_error > tolerance."""

    @given(
        pairs=layer_value_dicts(),
        tolerance=positive_tolerance,
    )
    @settings(max_examples=100)
    def test_flagged_iff_max_exceeds_tolerance(self, pairs, tolerance):
        """flagged is True if and only if max_relative_error > tolerance."""
        fort_dict, cpp_dict = pairs

        result = compute_error_metrics(fort_dict, cpp_dict, tolerance)

        assert result.flagged == (result.max_relative_error > tolerance)


class TestErrorMetricsPropertyZeroReferenceCount:
    """Property: num_zero_reference counts elements where fortran value was 0."""

    @given(
        pairs=layer_value_dicts(),
        tolerance=positive_tolerance,
    )
    @settings(max_examples=100)
    def test_num_zero_reference_count(self, pairs, tolerance):
        """num_zero_reference equals the count of elements with fort_val == 0."""
        fort_dict, cpp_dict = pairs

        result = compute_error_metrics(fort_dict, cpp_dict, tolerance)

        # Manually count zero reference values
        expected_count = 0
        common_layers = sorted(set(fort_dict.keys()) & set(cpp_dict.keys()))
        for layer in common_layers:
            fort_vals = fort_dict[layer]
            cpp_vals = cpp_dict[layer]
            num_species = min(len(fort_vals), len(cpp_vals))
            for i in range(num_species):
                if fort_vals[i] == 0.0:
                    expected_count += 1

        assert result.num_zero_reference == expected_count
