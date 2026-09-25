"""Profile condition management and scenario matrix generation.

Provides default atmospheric profile definitions, filtering helpers, and
the prepare_profile function that writes modified atmos_PTClds.dat files
by scaling cloud/aerosol fields.
"""

import os
import re
from typing import List, Optional

from benchmark.types import ProfileCondition, Scenario

# ---------------------------------------------------------------------------
# Default profile conditions and SZA values
# ---------------------------------------------------------------------------

DEFAULT_PROFILES: List[ProfileCondition] = [
    ProfileCondition(
        name="clear-sky",
        description="No clouds or aerosols",
        cloud_scale=0.0,
        aerosol_scale=0.0,
    ),
    ProfileCondition(
        name="cloudy",
        description="Standard cloud loading",
        cloud_scale=1.0,
        aerosol_scale=1.0,
    ),
    ProfileCondition(
        name="aerosol-loaded",
        description="Enhanced aerosol",
        cloud_scale=0.5,
        aerosol_scale=3.0,
    ),
]

DEFAULT_SZAS: List[float] = [0.0, 30.0, 60.0]

# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------


def get_profiles(subset: Optional[List[str]] = None) -> List[ProfileCondition]:
    """Return DEFAULT_PROFILES or a filtered subset matching given names.

    Args:
        subset: Optional list of profile names to include. If None, all
            default profiles are returned.

    Returns:
        List of ProfileCondition objects matching the requested names.

    Raises:
        ValueError: If a requested name does not match any default profile.
    """
    if subset is None:
        return list(DEFAULT_PROFILES)

    available = {p.name: p for p in DEFAULT_PROFILES}
    result: List[ProfileCondition] = []
    for name in subset:
        if name not in available:
            valid_names = ", ".join(sorted(available.keys()))
            raise ValueError(
                f"Unknown profile '{name}'. Available profiles: {valid_names}"
            )
        result.append(available[name])
    return result


def get_szas(subset: Optional[List[float]] = None) -> List[float]:
    """Return DEFAULT_SZAS or a filtered subset matching given values.

    Args:
        subset: Optional list of SZA values to include. If None, all
            default SZA values are returned.

    Returns:
        List of SZA float values matching the requested values.

    Raises:
        ValueError: If a requested SZA value is not in the default list.
    """
    if subset is None:
        return list(DEFAULT_SZAS)

    valid_set = set(DEFAULT_SZAS)
    result: List[float] = []
    for sza in subset:
        if sza not in valid_set:
            valid_vals = ", ".join(str(v) for v in DEFAULT_SZAS)
            raise ValueError(
                f"Unknown SZA value {sza}. Available values: {valid_vals}"
            )
        result.append(sza)
    return result


# ---------------------------------------------------------------------------
# Scenario matrix generation
# ---------------------------------------------------------------------------


def build_scenario_matrix(
    profiles: List[ProfileCondition],
    szas: List[float],
) -> List[Scenario]:
    """Build the Cartesian product of profiles and SZA values.

    Args:
        profiles: List of profile conditions.
        szas: List of solar zenith angle values in degrees.

    Returns:
        List of Scenario objects, one per (profile, sza) combination.
    """
    return [Scenario(profile=p, sza=s) for p in profiles for s in szas]


# ---------------------------------------------------------------------------
# Profile file preparation
# ---------------------------------------------------------------------------

# Pattern matching the atmospheric layer data lines (layer number + columns
# including AER-P values).  Example line:
#   1   0.0000000  1.0000000  299.5 0.80      21.84  0.000   0  0.000   0
_ATMOS_LAYER_RE = re.compile(
    r"^(\s*\d+)"  # layer number
    r"(\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+)"  # eta-A, eta-B, Temp, RH, Z(m)
    r"(\s+)(\S+)(\s+\S+\s+)(\S+)(\s+\S+)"  # AER-P1 NDA1 AER-P2 NDA2
    r"\s*$"
)

# Pattern matching cloud section data lines.  Example line:
#  29   0.40000E+00                               2.00000E-05   2.00000E-05
_CLOUD_LINE_RE = re.compile(
    r"^(\s*\d+)"  # layer number
    r"(\s+)(\S+)"  # cloud fraction
    r"(\s+)(\S+)"  # WLC
    r"(\s+)(\S+)"  # WIC
    r"\s*$"
)


def prepare_profile(profile: ProfileCondition, tables_dir: str) -> str:
    """Generate a modified atmos_PTClds.dat applying profile scaling.

    Reads the base atmos_PTClds.dat from *tables_dir*, applies the profile's
    cloud_scale and aerosol_scale multipliers to the relevant columns, and
    writes the modified file back to the same location.

    The aerosol columns (AER-P) in the atmospheric layer section are scaled
    by aerosol_scale.  The cloud section columns (cloud fraction, WLC, WIC)
    are scaled by cloud_scale.

    Args:
        profile: The profile condition defining scale factors.
        tables_dir: Path to the directory containing atmos_PTClds.dat.

    Returns:
        Path to the written (modified) atmos_PTClds.dat file.

    Raises:
        FileNotFoundError: If the base atmos_PTClds.dat does not exist.
    """
    filepath = os.path.join(tables_dir, "atmos_PTClds.dat")

    if not os.path.isfile(filepath):
        raise FileNotFoundError(
            f"Base profile file not found: {filepath}"
        )

    with open(filepath, "r") as f:
        lines = f.readlines()

    modified_lines: List[str] = []
    in_cloud_section = False

    for line in lines:
        # Detect start of cloud section
        if line.startswith("I=") and "Cloud" in line:
            in_cloud_section = True
            modified_lines.append(line)
            continue

        if in_cloud_section:
            m = _CLOUD_LINE_RE.match(line)
            if m:
                layer = m.group(1)
                sp1 = m.group(2)
                cloud_frac = m.group(3)
                sp2 = m.group(4)
                wlc = m.group(5)
                sp3 = m.group(6)
                wic = m.group(7)

                # Scale cloud fields
                new_cf = _scale_value(cloud_frac, profile.cloud_scale)
                new_wlc = _scale_value(wlc, profile.cloud_scale)
                new_wic = _scale_value(wic, profile.cloud_scale)

                modified_lines.append(
                    f"{layer}{sp1}{new_cf}{sp2}{new_wlc}{sp3}{new_wic}\n"
                )
            else:
                modified_lines.append(line)
        else:
            m = _ATMOS_LAYER_RE.match(line)
            if m:
                layer = m.group(1)
                mid_cols = m.group(2)
                sp_a1 = m.group(3)
                aer1 = m.group(4)
                mid_nda = m.group(5)
                aer2 = m.group(6)
                end = m.group(7)

                # Scale aerosol fields
                new_aer1 = _scale_float(aer1, profile.aerosol_scale)
                new_aer2 = _scale_float(aer2, profile.aerosol_scale)

                modified_lines.append(
                    f"{layer}{mid_cols}{sp_a1}{new_aer1}{mid_nda}{new_aer2}{end}\n"
                )
            else:
                modified_lines.append(line)

    with open(filepath, "w") as f:
        f.writelines(modified_lines)

    return filepath


def _scale_value(value_str: str, scale: float) -> str:
    """Scale a Fortran-format floating point value string.

    Preserves the original format (e.g., '0.40000E+00' stays in E notation).

    Args:
        value_str: Original value string (e.g., '2.00000E-05' or '0.40000E+00').
        scale: Multiplier to apply.

    Returns:
        Scaled value formatted to match the original width and style.
    """
    original_val = float(value_str)
    new_val = original_val * scale

    # Preserve Fortran scientific notation format
    if "E" in value_str.upper():
        # Match the original format width
        width = len(value_str)
        # Determine decimal places from original
        # e.g., '2.00000E-05' has 5 decimal places before E
        parts = value_str.upper().split("E")
        decimal_part = parts[0].split(".")
        if len(decimal_part) > 1:
            decimal_places = len(decimal_part[1])
        else:
            decimal_places = 5
        formatted = f"{new_val:.{decimal_places}E}"
        # Pad to match original width
        formatted = formatted.rjust(width)
        return formatted

    # Plain decimal format (e.g., '0.000')
    width = len(value_str)
    if "." in value_str:
        decimal_places = len(value_str.split(".")[1])
        formatted = f"{new_val:.{decimal_places}f}"
    else:
        formatted = f"{new_val:.3f}"
    return formatted.rjust(width)


def _scale_float(value_str: str, scale: float) -> str:
    """Scale a plain decimal float value and preserve width.

    Args:
        value_str: Original value string (e.g., '0.000').
        scale: Multiplier to apply.

    Returns:
        Scaled value formatted to match the original width.
    """
    original_val = float(value_str)
    new_val = original_val * scale
    width = len(value_str)
    if "." in value_str:
        decimal_places = len(value_str.split(".")[1])
        formatted = f"{new_val:.{decimal_places}f}"
    else:
        formatted = f"{new_val:.3f}"
    return formatted.rjust(width)
