"""Property-based tests for benchmark.report module.

**Validates: Requirements 6.2**

Property 6: Structured output field completeness
For any list of ScenarioResult objects, the JSON output produced by
write_json_output is valid JSON, contains a scenarios array with one entry
per result, and each entry contains all required fields.
"""

import json
import os
import tempfile

from hypothesis import given, settings
from hypothesis import strategies as st

from benchmark.report import write_json_output
from benchmark.types import ScenarioResult


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

# Strategy for generating ScenarioResult objects
scenario_result_st = st.builds(
    ScenarioResult,
    profile_name=st.text(
        min_size=1, max_size=20,
        alphabet=st.characters(whitelist_categories=("L", "N", "Pd")),
    ),
    sza=st.floats(min_value=0.0, max_value=90.0, allow_nan=False, allow_infinity=False),
    fortran_time=st.floats(min_value=0.0, max_value=1000.0, allow_nan=False, allow_infinity=False),
    cpp_time=st.floats(min_value=0.0, max_value=1000.0, allow_nan=False, allow_infinity=False),
    speedup=st.floats(min_value=0.0, max_value=1000.0, allow_nan=False, allow_infinity=False),
    max_relative_error=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
    mean_relative_error=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
    flagged=st.booleans(),
    fortran_failed=st.booleans(),
    cpp_failed=st.booleans(),
)

# Strategy for lists of ScenarioResult objects (including empty)
scenario_results_st = st.lists(scenario_result_st, min_size=0, max_size=10)

# Strategy for a positive tolerance
positive_tolerance = st.floats(
    min_value=1e-15, max_value=1e10, allow_nan=False, allow_infinity=False
)


# ---------------------------------------------------------------------------
# Property 6: Structured output field completeness
# ---------------------------------------------------------------------------


class TestJsonOutputFieldCompleteness:
    """**Validates: Requirements 6.2**"""

    @given(results=scenario_results_st, tolerance=positive_tolerance)
    @settings(max_examples=200)
    def test_output_is_valid_json(self, results, tolerance):
        """write_json_output produces a file containing valid JSON."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            output_path = f.name

        try:
            write_json_output(results, output_path, tolerance)

            with open(output_path, "r") as f:
                data = json.load(f)

            assert isinstance(data, dict)
        finally:
            os.unlink(output_path)

    @given(results=scenario_results_st, tolerance=positive_tolerance)
    @settings(max_examples=200)
    def test_scenarios_array_length_matches_results(self, results, tolerance):
        """The JSON scenarios array has one entry per ScenarioResult."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            output_path = f.name

        try:
            write_json_output(results, output_path, tolerance)

            with open(output_path, "r") as f:
                data = json.load(f)

            assert "scenarios" in data
            assert isinstance(data["scenarios"], list)
            assert len(data["scenarios"]) == len(results)
        finally:
            os.unlink(output_path)

    @given(results=scenario_results_st, tolerance=positive_tolerance)
    @settings(max_examples=200)
    def test_each_entry_contains_all_required_fields(self, results, tolerance):
        """Each scenario entry contains all required fields."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            output_path = f.name

        try:
            write_json_output(results, output_path, tolerance)

            with open(output_path, "r") as f:
                data = json.load(f)

            required_fields = {
                "profile",
                "sza",
                "fortran_elapsed_s",
                "cpp_elapsed_s",
                "speedup",
                "max_relative_error",
                "mean_relative_error",
                "status",
            }

            for entry in data["scenarios"]:
                assert required_fields.issubset(entry.keys()), (
                    f"Missing fields: {required_fields - set(entry.keys())}"
                )
        finally:
            os.unlink(output_path)

    @given(results=scenario_results_st, tolerance=positive_tolerance)
    @settings(max_examples=200)
    def test_summary_object_contains_required_fields(self, results, tolerance):
        """The JSON summary object contains all required summary fields."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            output_path = f.name

        try:
            write_json_output(results, output_path, tolerance)

            with open(output_path, "r") as f:
                data = json.load(f)

            assert "summary" in data
            summary = data["summary"]

            required_summary_fields = {
                "overall_max_error",
                "overall_mean_speedup",
                "total_scenarios",
                "flagged_scenarios",
            }

            assert required_summary_fields.issubset(summary.keys()), (
                f"Missing summary fields: {required_summary_fields - set(summary.keys())}"
            )
        finally:
            os.unlink(output_path)

    @given(results=scenario_results_st, tolerance=positive_tolerance)
    @settings(max_examples=200)
    def test_status_field_reflects_flagged(self, results, tolerance):
        """The status field is 'pass' when not flagged and 'fail' when flagged."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            output_path = f.name

        try:
            write_json_output(results, output_path, tolerance)

            with open(output_path, "r") as f:
                data = json.load(f)

            for i, entry in enumerate(data["scenarios"]):
                expected_status = "fail" if results[i].flagged else "pass"
                assert entry["status"] == expected_status, (
                    f"Scenario {i}: expected status '{expected_status}', "
                    f"got '{entry['status']}' (flagged={results[i].flagged})"
                )
        finally:
            os.unlink(output_path)
