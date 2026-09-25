#!/usr/bin/env python3
"""Cloud-J Fortran vs C++ Benchmark Harness.

Main orchestration script that runs both Fortran and C++ Cloud-J standalone
executables across a scenario matrix of atmospheric profiles and solar zenith
angles, compares their J-value outputs for numerical accuracy, measures
wall-clock performance, and generates terminal and structured reports.

Usage:
    python -m benchmark.run_benchmark \
        --fortran bin/cloudj_standalone \
        --cpp bin/cloudj_standalone_cpp \
        --profiles clear-sky cloudy \
        --tolerance 1e-6 \
        --output results/benchmark.json
"""

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import time
from typing import List, Optional

from benchmark.metrics import compute_error_metrics, compute_speedup
from benchmark.parser import parse_jvalues
from benchmark.profiles import build_scenario_matrix, get_profiles, get_szas, prepare_profile
from benchmark.report import format_terminal_report, write_csv_output, write_json_output
from benchmark.types import ExecutionResult, ScenarioResult


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the benchmark harness.

    Returns:
        Namespace with all parsed CLI options.
    """
    parser = argparse.ArgumentParser(
        description="Cloud-J Fortran vs C++ Benchmark"
    )
    parser.add_argument(
        "--fortran", required=True, help="Path to cloudj_standalone"
    )
    parser.add_argument(
        "--cpp", required=True, help="Path to cloudj_standalone_cpp"
    )
    parser.add_argument(
        "--profiles",
        nargs="*",
        default=None,
        help="Subset of profile conditions to run",
    )
    parser.add_argument(
        "--sza",
        nargs="*",
        type=float,
        default=None,
        help="Subset of SZA values to run",
    )
    parser.add_argument(
        "--output",
        default="benchmark_results.json",
        help="Output file path for structured results",
    )
    parser.add_argument(
        "--format",
        choices=["json", "csv"],
        default="json",
        help="Output format",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="Max relative error threshold for flagging divergence",
    )
    parser.add_argument(
        "--benchmark-iters",
        action="store_true",
        help="Use --benchmark flag on C++ executable for throughput timing",
    )
    return parser.parse_args()


def run_executable(
    exe_path: str,
    working_dir: str,
    extra_args: Optional[List[str]] = None,
) -> ExecutionResult:
    """Run an executable, capture output, and measure wall-clock time.

    Uses a shell wrapper with 'ulimit -s hard' to prevent Fortran
    stack overflow segfaults on macOS (where 'unlimited' is rejected).
    """
    cmd_parts = [exe_path] + (extra_args or [])
    # Wrap in shell with maximum allowed stack size for Fortran compatibility.
    # On macOS, 'ulimit -s unlimited' fails even with 2>/dev/null suppression,
    # so we use 'ulimit -s hard' which sets the soft limit to the hard limit.
    shell_cmd = "ulimit -s hard 2>/dev/null; exec " + " ".join(
        shlex.quote(p) for p in cmd_parts
    )

    start = time.perf_counter()
    result = subprocess.run(
        ["bash", "-c", shell_cmd],
        cwd=working_dir,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - start
    return ExecutionResult(
        stdout=result.stdout,
        elapsed_seconds=elapsed,
        returncode=result.returncode,
        error_message=result.stderr if result.returncode != 0 else None,
    )


def main() -> None:
    """Main entry point for the benchmark harness."""
    args = parse_args()

    # Validate executable paths exist
    fortran_path = os.path.abspath(args.fortran)
    cpp_path = os.path.abspath(args.cpp)

    if not os.path.isfile(fortran_path):
        print(
            f"Error: Fortran executable not found: {fortran_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.isfile(cpp_path):
        print(
            f"Error: C++ executable not found: {cpp_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Determine the working directory: the directory containing the executables
    # and the tables/ directory. Use the Fortran executable's directory.
    working_dir = os.path.dirname(fortran_path)

    # Determine tables directory for profile preparation
    tables_dir = os.path.join(working_dir, "tables")

    # Get profiles and SZA values (full or subset)
    profiles = get_profiles(args.profiles)
    szas = get_szas(args.sza)

    # Build scenario matrix
    scenarios = build_scenario_matrix(profiles, szas)

    # Determine number of SZA blocks per run (executables loop over all SZAs)
    num_sza_blocks = len(szas)

    # Collect results
    results: List[ScenarioResult] = []

    # Group scenarios by profile (each executable run covers all SZAs for a profile)
    seen_profiles = []
    for profile in profiles:
        if profile.name not in [p.name for p in seen_profiles]:
            seen_profiles.append(profile)

    # Back up the original atmos_PTClds.dat before any profile preparation
    # mutates it in place. This file is shared with the standalone executables
    # and the ctest comparison tests, so we must always restore it afterward.
    profile_path = os.path.join(tables_dir, "atmos_PTClds.dat")
    backup_path = None
    if os.path.isfile(profile_path):
        backup_path = profile_path + ".benchmark_backup"
        shutil.copy2(profile_path, backup_path)

    try:
        for profile in seen_profiles:
            # Prepare the profile file (modify atmos_PTClds.dat).  Restore the
            # pristine backup first: prepare_profile scales the *current* file
            # in place, so without this the previous profile's scaling (e.g.
            # clear-sky zeroing the cloud columns) would compound into the
            # next one and mask the cloudy scenarios.
            if backup_path is not None and os.path.isfile(backup_path):
                shutil.copy2(backup_path, profile_path)
            try:
                prepare_profile(profile, tables_dir)
            except FileNotFoundError as e:
                print(f"Warning: Could not prepare profile '{profile.name}': {e}", file=sys.stderr)
                # Add failed results for all SZAs in this profile
                for sza in szas:
                    results.append(ScenarioResult(
                        profile_name=profile.name,
                        sza=sza,
                        fortran_time=0.0,
                        cpp_time=0.0,
                        speedup=0.0,
                        max_relative_error=0.0,
                        mean_relative_error=0.0,
                        flagged=False,
                        fortran_failed=True,
                        cpp_failed=True,
                    ))
                continue

            # Run Fortran executable
            fortran_result = run_executable(fortran_path, working_dir)

            # Run C++ executable (with --benchmark flag if requested)
            cpp_extra_args = ["--benchmark"] if args.benchmark_iters else None
            cpp_result = run_executable(cpp_path, working_dir, extra_args=cpp_extra_args)

            # Handle non-zero exit codes
            fortran_failed = fortran_result.returncode != 0
            cpp_failed = cpp_result.returncode != 0

            if fortran_failed:
                print(
                    f"Warning: Fortran executable failed for profile '{profile.name}' "
                    f"(exit code {fortran_result.returncode}): "
                    f"{fortran_result.error_message}",
                    file=sys.stderr,
                )

            if cpp_failed:
                print(
                    f"Warning: C++ executable failed for profile '{profile.name}' "
                    f"(exit code {cpp_result.returncode}): "
                    f"{cpp_result.error_message}",
                    file=sys.stderr,
                )

            # If both failed, record failure for all SZAs and continue
            if fortran_failed and cpp_failed:
                for sza in szas:
                    results.append(ScenarioResult(
                        profile_name=profile.name,
                        sza=sza,
                        fortran_time=0.0,
                        cpp_time=0.0,
                        speedup=0.0,
                        max_relative_error=0.0,
                        mean_relative_error=0.0,
                        flagged=False,
                        fortran_failed=True,
                        cpp_failed=True,
                    ))
                continue

            # Parse J-values from outputs
            fortran_blocks = parse_jvalues(fortran_result.stdout) if not fortran_failed else []
            cpp_blocks = parse_jvalues(cpp_result.stdout) if not cpp_failed else []

            # Compute per-SZA timing (divide total time by number of SZA blocks)
            fortran_per_sza = (
                fortran_result.elapsed_seconds / num_sza_blocks
                if num_sza_blocks > 0 and not fortran_failed
                else 0.0
            )
            cpp_per_sza = (
                cpp_result.elapsed_seconds / num_sza_blocks
                if num_sza_blocks > 0 and not cpp_failed
                else 0.0
            )

            # Process each SZA block
            for sza_idx, sza in enumerate(szas):
                # Get the corresponding J-value blocks
                fortran_jvals = (
                    fortran_blocks[sza_idx].values
                    if sza_idx < len(fortran_blocks)
                    else {}
                )
                cpp_jvals = (
                    cpp_blocks[sza_idx].values
                    if sza_idx < len(cpp_blocks)
                    else {}
                )

                # Compute error metrics if we have data from both
                if fortran_jvals and cpp_jvals and not fortran_failed and not cpp_failed:
                    error_metrics = compute_error_metrics(
                        fortran_jvals, cpp_jvals, args.tolerance
                    )
                    speedup = compute_speedup(fortran_per_sza, cpp_per_sza)
                else:
                    from benchmark.types import ErrorMetrics

                    error_metrics = ErrorMetrics(
                        max_relative_error=0.0,
                        mean_relative_error=0.0,
                        num_elements=0,
                        num_zero_reference=0,
                        flagged=False,
                    )
                    speedup = 0.0

                results.append(ScenarioResult(
                    profile_name=profile.name,
                    sza=sza,
                    fortran_time=fortran_per_sza,
                    cpp_time=cpp_per_sza,
                    speedup=speedup,
                    max_relative_error=error_metrics.max_relative_error,
                    mean_relative_error=error_metrics.mean_relative_error,
                    flagged=error_metrics.flagged,
                    fortran_failed=fortran_failed,
                    cpp_failed=cpp_failed,
                ))

        # Generate terminal report
        report = format_terminal_report(results, args.tolerance)
        print(report)

        # Write structured output
        if args.format == "json":
            write_json_output(results, args.output, tolerance=args.tolerance)
        else:
            write_csv_output(results, args.output)

        print(f"\nResults written to: {args.output}")
    finally:
        # Always restore the original atmos_PTClds.dat so we never leave the
        # shared table file mutated (it is also read by the standalone
        # executables and the ctest comparison tests).
        if backup_path is not None and os.path.isfile(backup_path):
            shutil.move(backup_path, profile_path)


if __name__ == "__main__":
    main()
