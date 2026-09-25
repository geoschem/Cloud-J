"""Report generation for the Cloud-J benchmark harness.

Provides functions to format benchmark results as a terminal table,
write structured JSON output, and write CSV output for CI consumption.
"""

import csv
import json
import os
from datetime import datetime, timezone
from typing import List

from benchmark.types import ScenarioResult


def format_terminal_report(results: List[ScenarioResult], tolerance: float) -> str:
    """Format results as a human-readable fixed-width terminal table.

    Includes one row per scenario with all fields, plus a summary line
    with overall max error, mean speedup, and flagged count.
    Flagged rows are prefixed with a warning indicator (⚠).

    Args:
        results: List of ScenarioResult objects from the benchmark run.
        tolerance: The tolerance threshold used for flagging.

    Returns:
        Formatted multi-line string suitable for printing to stdout.
    """
    lines: List[str] = []

    # Header
    tol_str = f"{tolerance:.1e}"
    lines.append(f"Cloud-J Benchmark Results (tolerance: {tol_str})")
    lines.append("\u2550" * 79)

    # Column headers
    header = (
        f"{'Profile':<16} {'SZA':>5}  "
        f"{'Fortran(s)':>10}  {'C++(s)':>8}  "
        f"{'Speedup':>8}  {'MaxRelErr':>11}  {'MeanRelErr':>11}"
    )
    lines.append(header)
    lines.append("\u2500" * 79)

    # Data rows
    for r in results:
        prefix = "\u26a0 " if r.flagged else "  "
        profile_field = f"{prefix}{r.profile_name}"

        row = (
            f"{profile_field:<18} {r.sza:>5.1f}  "
            f"{r.fortran_time:>10.3f}  {r.cpp_time:>8.3f}  "
            f"{r.speedup:>7.2f}x  "
            f"{r.max_relative_error:>11.2e}  {r.mean_relative_error:>11.2e}"
        )
        lines.append(row)

    # Summary separator
    lines.append("\u2500" * 79)

    # Compute summary statistics
    if results:
        overall_max_error = max(r.max_relative_error for r in results)
        overall_mean_speedup = (
            sum(r.speedup for r in results) / len(results)
        )
        flagged_count = sum(1 for r in results if r.flagged)
        total_count = len(results)
    else:
        overall_max_error = 0.0
        overall_mean_speedup = 0.0
        flagged_count = 0
        total_count = 0

    summary = (
        f"SUMMARY: Max Error = {overall_max_error:.2e} | "
        f"Mean Speedup = {overall_mean_speedup:.2f}x | "
        f"Flagged: {flagged_count}/{total_count}"
    )
    lines.append(summary)
    lines.append("\u2550" * 79)

    return "\n".join(lines)


def write_json_output(
    results: List[ScenarioResult], output_path: str, tolerance: float = 1e-6
) -> None:
    """Write structured benchmark results to a JSON file.

    Creates parent directories if they don't exist. The JSON includes
    a timestamp, tolerance, per-scenario data, and an overall summary.

    Args:
        results: List of ScenarioResult objects from the benchmark run.
        output_path: File path for the JSON output.
        tolerance: The tolerance threshold used for flagging.
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Build scenarios array
    scenarios = []
    for r in results:
        scenarios.append({
            "profile": r.profile_name,
            "sza": r.sza,
            "fortran_elapsed_s": r.fortran_time,
            "cpp_elapsed_s": r.cpp_time,
            "speedup": r.speedup,
            "max_relative_error": r.max_relative_error,
            "mean_relative_error": r.mean_relative_error,
            "status": "fail" if r.flagged else "pass",
        })

    # Compute summary
    if results:
        overall_max_error = max(r.max_relative_error for r in results)
        overall_mean_speedup = sum(r.speedup for r in results) / len(results)
        flagged_count = sum(1 for r in results if r.flagged)
    else:
        overall_max_error = 0.0
        overall_mean_speedup = 0.0
        flagged_count = 0

    output = {
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tolerance": tolerance,
        "scenarios": scenarios,
        "summary": {
            "overall_max_error": overall_max_error,
            "overall_mean_speedup": overall_mean_speedup,
            "total_scenarios": len(results),
            "flagged_scenarios": flagged_count,
        },
    }

    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)


def write_csv_output(results: List[ScenarioResult], output_path: str) -> None:
    """Write structured benchmark results to a CSV file.

    Creates parent directories if they don't exist.

    Headers: profile,sza,fortran_elapsed_s,cpp_elapsed_s,speedup,
             max_relative_error,mean_relative_error,status

    Args:
        results: List of ScenarioResult objects from the benchmark run.
        output_path: File path for the CSV output.
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    fieldnames = [
        "profile",
        "sza",
        "fortran_elapsed_s",
        "cpp_elapsed_s",
        "speedup",
        "max_relative_error",
        "mean_relative_error",
        "status",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            writer.writerow({
                "profile": r.profile_name,
                "sza": r.sza,
                "fortran_elapsed_s": r.fortran_time,
                "cpp_elapsed_s": r.cpp_time,
                "speedup": r.speedup,
                "max_relative_error": r.max_relative_error,
                "mean_relative_error": r.mean_relative_error,
                "status": "fail" if r.flagged else "pass",
            })
