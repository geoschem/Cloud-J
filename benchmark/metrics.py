"""Error and speedup computation for the Cloud-J benchmark harness.

Provides functions to compute relative/absolute error metrics between
Fortran and C++ J-value outputs, and speedup ratios from timing data.
"""

from typing import Dict, List

from benchmark.types import ErrorMetrics


def compute_error_metrics(
    fortran_values: Dict[int, List[float]],
    cpp_values: Dict[int, List[float]],
    tolerance: float,
) -> ErrorMetrics:
    """Compute relative error metrics between Fortran and C++ J-values.

    For each (layer, species) pair present in both dictionaries:
    - If fortran_value == 0: use absolute difference |cpp - fortran|
    - Otherwise: use relative error |cpp - fortran| / |fortran|

    Args:
        fortran_values: Mapping of layer_number to list of J-values per species
            from the Fortran executable output.
        cpp_values: Mapping of layer_number to list of J-values per species
            from the C++ executable output.
        tolerance: Maximum relative error threshold for flagging divergence.

    Returns:
        ErrorMetrics with max/mean error, element counts, and flagging status.
    """
    errors: List[float] = []
    num_zero_reference = 0

    # Compare only layers present in both outputs
    common_layers = sorted(set(fortran_values.keys()) & set(cpp_values.keys()))

    for layer in common_layers:
        fort_vals = fortran_values[layer]
        cpp_vals = cpp_values[layer]

        # Compare only the species columns present in both
        num_species = min(len(fort_vals), len(cpp_vals))

        for i in range(num_species):
            fort_val = fort_vals[i]
            cpp_val = cpp_vals[i]

            if fort_val == 0.0:
                # Use absolute difference when reference is zero
                error = abs(cpp_val - fort_val)
                num_zero_reference += 1
            else:
                # Use relative error
                error = abs(cpp_val - fort_val) / abs(fort_val)

            errors.append(error)

    num_elements = len(errors)

    if num_elements == 0:
        return ErrorMetrics(
            max_relative_error=0.0,
            mean_relative_error=0.0,
            num_elements=0,
            num_zero_reference=0,
            flagged=False,
        )

    max_error = max(errors)
    mean_error = sum(errors) / num_elements

    return ErrorMetrics(
        max_relative_error=max_error,
        mean_relative_error=mean_error,
        num_elements=num_elements,
        num_zero_reference=num_zero_reference,
        flagged=max_error > tolerance,
    )


def compute_speedup(fortran_time: float, cpp_time: float) -> float:
    """Compute speedup ratio: fortran_time / cpp_time.

    Args:
        fortran_time: Elapsed time in seconds for the Fortran executable.
        cpp_time: Elapsed time in seconds for the C++ executable.

    Returns:
        Speedup ratio. Returns float('inf') if cpp_time is zero or negative.
    """
    if cpp_time <= 0.0:
        return float('inf')
    return fortran_time / cpp_time
