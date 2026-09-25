"""Shared data types for the Cloud-J benchmark harness.

All dataclasses used across the benchmark modules are defined here to
avoid circular imports and provide a single source of truth for the
data model.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class ProfileCondition:
    """An atmospheric column scenario defined by cloud/aerosol scaling.

    Attributes:
        name: Short identifier, e.g. "clear-sky", "cloudy", "aerosol-loaded".
        description: Human-readable description of the condition.
        cloud_scale: Multiplier for cloud fields (0.0 = clear sky).
        aerosol_scale: Multiplier for aerosol fields.
    """

    name: str
    description: str
    cloud_scale: float
    aerosol_scale: float


@dataclass
class Scenario:
    """A single benchmark test case combining a profile and SZA.

    Attributes:
        profile: The atmospheric profile condition for this scenario.
        sza: Solar Zenith Angle in degrees.
    """

    profile: ProfileCondition
    sza: float


@dataclass
class ExecutionResult:
    """Result of running a single executable.

    Attributes:
        stdout: Captured standard output from the process.
        elapsed_seconds: Wall-clock time in seconds.
        returncode: Process exit code (0 = success).
        error_message: Captured stderr if the process failed, else None.
    """

    stdout: str
    elapsed_seconds: float
    returncode: int
    error_message: Optional[str] = None


@dataclass
class JValueBlock:
    """Parsed J-value output for a single SZA block.

    Attributes:
        sza_index: 0-based index of this SZA block in the output.
        species: Column headers (species names) from the output table.
        values: Mapping of layer_number to list of J-values per species.
    """

    sza_index: int
    species: List[str] = field(default_factory=list)
    values: Dict[int, List[float]] = field(default_factory=dict)


@dataclass
class ErrorMetrics:
    """Numerical accuracy metrics comparing Fortran and C++ J-values.

    Attributes:
        max_relative_error: Maximum relative (or absolute) error across all elements.
        mean_relative_error: Mean relative (or absolute) error across all elements.
        num_elements: Total number of elements compared.
        num_zero_reference: Count of elements where absolute diff was used
            (because the Fortran reference value was zero).
        flagged: True if max_relative_error exceeds the tolerance threshold.
    """

    max_relative_error: float
    mean_relative_error: float
    num_elements: int
    num_zero_reference: int
    flagged: bool


@dataclass
class ScenarioResult:
    """Aggregated result for a single benchmark scenario.

    Attributes:
        profile_name: Name of the profile condition.
        sza: Solar Zenith Angle in degrees.
        fortran_time: Fortran executable elapsed time in seconds.
        cpp_time: C++ executable elapsed time in seconds.
        speedup: Ratio of fortran_time / cpp_time.
        max_relative_error: Maximum relative error for this scenario.
        mean_relative_error: Mean relative error for this scenario.
        flagged: True if max error exceeds tolerance.
        fortran_failed: True if the Fortran executable returned non-zero.
        cpp_failed: True if the C++ executable returned non-zero.
    """

    profile_name: str
    sza: float
    fortran_time: float
    cpp_time: float
    speedup: float
    max_relative_error: float
    mean_relative_error: float
    flagged: bool
    fortran_failed: bool
    cpp_failed: bool
