# Cloud-J Benchmark Harness

Compares the C++ Cloud-J port (`cloudj_standalone_cpp`) against the Fortran reference implementation (`cloudj_standalone`) across a matrix of atmospheric profiles and solar zenith angles. Reports wall-clock performance and numerical accuracy.

## Prerequisites

1. **Both executables must be built:**
   - `cloudj_standalone` (Fortran) — built from the standard Fortran source
   - `cloudj_standalone_cpp` (C++) — built with `CLOUDJ_CPP_PORT=ON`

   Both are produced by the CMake build system when the C++ port is enabled, and land in `bin/` (the CMake `RUNTIME_OUTPUT_DIRECTORY` for this project) alongside a `bin/tables/` directory copied from the project's `tables/`.

2. **The `tables/` directory must be accessible** from the working directory where the executables run — normally `bin/tables/`, next to both executables.

3. **Python 3.6+** — the harness uses only the standard library (no third-party dependencies).

4. **macOS note**: the Fortran executable can segfault under the default 8MB stack limit due to large stack-allocated arrays. The harness's `run_executable()` already works around this automatically by running both executables through a shell wrapper that does `ulimit -s hard` before exec (macOS rejects `ulimit -s unlimited` outright for the stack resource, unlike Linux, so the harness raises the soft limit to whatever the OS reports as the hard limit instead). You don't need to do anything for this — it's mentioned here so a segfault from running the Fortran binary *directly* (bypassing the harness) isn't surprising.

## Usage

```bash
python benchmark/run_benchmark.py \
    --fortran <path-to-cloudj_standalone> \
    --cpp <path-to-cloudj_standalone_cpp>
```

Or as a module:

```bash
python -m benchmark.run_benchmark \
    --fortran bin/cloudj_standalone \
    --cpp bin/cloudj_standalone_cpp
```

## CLI Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--fortran` | Yes | — | Path to the Fortran `cloudj_standalone` executable |
| `--cpp` | Yes | — | Path to the C++ `cloudj_standalone_cpp` executable |
| `--profiles` | No | all | Subset of profile conditions to run |
| `--sza` | No | all | Subset of SZA values to run |
| `--output` | No | `benchmark_results.json` | Output file path for structured results |
| `--format` | No | `json` | Output format: `json` or `csv` |
| `--tolerance` | No | `1e-6` | Relative error threshold for flagging divergence |
| `--benchmark-iters` | No | off | Passes a bare `--benchmark` flag to the C++ executable, which the standalone driver treats as an alias for `--bench-iters 1000` (see the standalone driver's own in-memory benchmark mode below) |

### Profile conditions

Available profiles: `clear-sky`, `cloudy`, `aerosol-loaded`

### SZA values

Available SZA values: `0`, `30`, `60` (degrees)

### A note on the C++ executable's own `--bench-iters` flag

Separately from this Python harness, `cloudj_standalone_cpp` has its own `--bench-iters N` flag that loops the corrected `engine.cloud_jx()` call N times in-process (excluding one-time setup), useful for measuring pure in-memory per-call latency without process-startup/table-loading overhead. This Python harness's `--benchmark-iters` flag is a convenience wrapper that triggers a small (1000-iteration) version of that same in-memory loop on the C++ side via the `--benchmark` alias; it does not affect the Fortran side and does not change what's compared for numerical accuracy. If you want a dedicated in-memory throughput number rather than the full wall-clock comparison this harness produces, run `bin/cloudj_standalone_cpp --bench-iters 20000` directly.

## Examples

### Full benchmark run (all profiles, all SZAs)

```bash
python -m benchmark.run_benchmark \
    --fortran bin/cloudj_standalone \
    --cpp bin/cloudj_standalone_cpp
```

### Subset run (specific profiles and SZAs)

```bash
python -m benchmark.run_benchmark \
    --fortran bin/cloudj_standalone \
    --cpp bin/cloudj_standalone_cpp \
    --profiles clear-sky cloudy \
    --sza 0 60
```

### CSV output with custom tolerance

```bash
python -m benchmark.run_benchmark \
    --fortran bin/cloudj_standalone \
    --cpp bin/cloudj_standalone_cpp \
    --format csv \
    --output results/benchmark.csv \
    --tolerance 1e-8
```

### C++ internal benchmark mode

```bash
python -m benchmark.run_benchmark \
    --fortran bin/cloudj_standalone \
    --cpp bin/cloudj_standalone_cpp \
    --benchmark-iters
```

## Output

### Terminal Report

The harness prints a summary table to stdout after each run. This is real output from a run against the current build on a development machine (not illustrative placeholder numbers):

```
Cloud-J Benchmark Results (tolerance: 1.0e-06)
═══════════════════════════════════════════════════════════════════════════════
Profile            SZA  Fortran(s)    C++(s)   Speedup    MaxRelErr   MeanRelErr
───────────────────────────────────────────────────────────────────────────────
  clear-sky          0.0       0.079     0.019     4.07x     0.00e+00     0.00e+00
  clear-sky         30.0       0.079     0.019     4.07x     0.00e+00     0.00e+00
  clear-sky         60.0       0.079     0.019     4.07x     0.00e+00     0.00e+00
  cloudy             0.0       0.017     0.016     1.06x     0.00e+00     0.00e+00
  cloudy            30.0       0.017     0.016     1.06x     0.00e+00     0.00e+00
  cloudy            60.0       0.017     0.016     1.06x     0.00e+00     0.00e+00
  aerosol-loaded     0.0       0.017     0.017     1.03x     0.00e+00     0.00e+00
  aerosol-loaded    30.0       0.017     0.017     1.03x     0.00e+00     0.00e+00
  aerosol-loaded    60.0       0.017     0.017     1.03x     0.00e+00     0.00e+00
───────────────────────────────────────────────────────────────────────────────
SUMMARY: Max Error = 0.00e+00 | Mean Speedup = 2.05x | Flagged: 0/9
═══════════════════════════════════════════════════════════════════════════════
```

**Note on the first row's speedup**: `clear-sky` at SZA 0.0 above shows an inflated 4.07x speedup because that profile ran first, and the Fortran binary's first invocation in a run is consistently slower (0.079s here vs the ~0.017s steady-state seen in every subsequent row) — most likely OS-level cold-start effects (page faults, filesystem cache) rather than anything about the Fortran code itself. **Steady-state speedup is consistently around 1.0-1.1x** (rough parity, C++ marginally faster) — treat any single run's summary "Mean Speedup" with this in mind, since one slow first row skews the average (2.05x here). Run the harness 2-3 times and look at the non-first rows if you want a representative number; see `docs/C++_Header_Port_Completion_Report.md` for a fuller discussion of steady-state vs cold-start timing and the in-memory `--bench-iters` numbers used to characterize actual per-call performance.

### JSON Output

The default structured output format. Contains a timestamp, tolerance, per-scenario data, and a summary — this is the real output produced by the run shown above:

```json
{
  "timestamp": "2026-08-06T14:54:42Z",
  "tolerance": 1e-06,
  "scenarios": [
    {"profile": "clear-sky", "sza": 0.0, "fortran_elapsed_s": 0.0790997776, "cpp_elapsed_s": 0.0194538053, "speedup": 4.0660311094, "max_relative_error": 0.0, "mean_relative_error": 0.0, "status": "pass"},
    {"profile": "cloudy", "sza": 0.0, "fortran_elapsed_s": 0.0169765693, "cpp_elapsed_s": 0.0160614307, "speedup": 1.0569774045, "max_relative_error": 0.0, "mean_relative_error": 0.0, "status": "pass"},
    {"profile": "aerosol-loaded", "sza": 0.0, "fortran_elapsed_s": 0.0171444723, "cpp_elapsed_s": 0.0167066527, "speedup": 1.0262063071, "max_relative_error": 0.0, "mean_relative_error": 0.0, "status": "pass"}
  ],
  "summary": {
    "overall_max_error": 0.0,
    "overall_mean_speedup": 2.0497382737,
    "total_scenarios": 9,
    "flagged_scenarios": 0
  }
}
```
(Truncated to one representative row per profile above for brevity — a real run produces all 9 scenario entries, one per profile × SZA combination.)

### CSV Output

Flat tabular format with one row per scenario:

```
profile,sza,fortran_elapsed_s,cpp_elapsed_s,speedup,max_relative_error,mean_relative_error,status
clear-sky,0.0,0.0790997776,0.0194538053,4.0660311094,0.0,0.0,pass
cloudy,0.0,0.0169765693,0.0160614307,1.0569774045,0.0,0.0,pass
aerosol-loaded,0.0,0.0171444723,0.0167066527,1.0262063071,0.0,0.0,pass
```

## Interpreting Results

### Speedup

The speedup ratio is `Fortran time / C++ time`. A value greater than 1.0 means the C++ implementation is faster than Fortran for that run. As noted above, the first scenario run in any harness invocation tends to show an inflated speedup from Fortran cold-start effects — steady-state (later rows) speedup for the current C++ port is roughly 1.0-1.1x, not a large multiplier. See `docs/C++_Header_Port_Completion_Report.md` for the full performance analysis, including in-memory (non-process-startup) throughput numbers from the C++ executable's own `--bench-iters` mode.

### Flagged Scenarios

A scenario is flagged (marked with ⚠ in the terminal report, `"status": "fail"` in JSON) when its **maximum relative error exceeds the tolerance threshold**. This indicates a numerical divergence between the Fortran and C++ implementations for that particular profile/SZA combination.

As of the current build, all 9 scenarios pass with exactly `0.00e+00` max and mean relative error — the C++ port and Fortran reference produce numerically identical J-values at the precision this harness's text-based comparison can measure (see the caveat in `docs/C++_Header_Port_Completion_Report.md` about `e9.2`-format text comparison vs raw double precision). If you see flagged scenarios after making a code change, they may indicate:
- A bug introduced in the C++ port
- Differences in floating-point evaluation order
- Uninitialized data or missing initialization in the C++ code
- A mismatch between on-disk override tables and embedded tables, if using the C++ engine's `--tables-dir` option (see `include/cloudj/init.hpp` and the standalone driver's `--tables-dir` flag)

### Tolerance

The `--tolerance` flag sets the relative error threshold (default: `1e-6`). The relative error for each J-value element is computed as:

- `|cpp_value - fortran_value| / |fortran_value|` when the Fortran value is non-zero
- `|cpp_value - fortran_value|` (absolute difference) when the Fortran value is zero

A tighter tolerance (e.g., `1e-10`) will flag more scenarios; a looser tolerance (e.g., `1e-4`) will be more permissive.

### Error Metrics

- **Max Relative Error**: The worst-case divergence across all layers and species in a scenario. Use this to identify the most problematic J-value.
- **Mean Relative Error**: The average divergence. A low mean with a high max suggests an isolated issue rather than systematic drift.
