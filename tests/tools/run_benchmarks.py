# tests/tools/run_benchmarks.py
import sys
import os
import subprocess
import time
import argparse

def execute_binary_repeatedly(binary_path, iterations, use_headers_flag=True):
    """
    Measures the high-resolution execution time of a standalone binary over N iterations.
    """
    args = ["./" + os.path.basename(binary_path)]
    if use_headers_flag and "standalone_cpp" in binary_path:
        args.append("--use-headers-only")

    start_time = time.perf_counter()
    for i in range(1, iterations + 1):
        res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.path.dirname(binary_path))
        # Ignore returncode for Fortran standalone on Darwin/macOS if it successfully ran and printed results
        is_fortran = "standalone_cpp" not in binary_path
        if res.returncode != 0 and not (is_fortran and sys.platform == "darwin"):
            print(f"Error: Binary {binary_path} failed to execute inside {os.path.dirname(binary_path)}.", file=sys.stderr)
            return None
        if i % 500 == 0:
            print(f"  -> Completed {i} / {iterations} runs...", flush=True)
    end_time = time.perf_counter()
    return end_time - start_time

def run_benchmarks(iterations):
    print("=================================================================")
    print(f"Starting High-Resolution Performance Benchmarking ({iterations} runs)...")
    print("=================================================================\n")

    fortran_bin = "./cloudj_standalone"
    cpp_bin = "./cloudj_standalone_cpp"

    # Locate paths dynamically
    if not os.path.exists(fortran_bin) or not os.path.exists(cpp_bin):
        for path in ["build/bin", "bin", "."]:
            f_path = os.path.join(path, "cloudj_standalone")
            c_path = os.path.join(path, "cloudj_standalone_cpp")
            if os.path.exists(f_path) and os.path.exists(c_path):
                fortran_bin = f_path
                cpp_bin = c_path
                break
        else:
            print("Error: Standalone binaries not found.", file=sys.stderr)
            return False

    # 1. Benchmark compiled Fortran version
    print("Benchmarking Fortran Standalone...")
    t_fortran = execute_binary_repeatedly(fortran_bin, iterations, use_headers_flag=False)
    if t_fortran is None: return False

    # 2. Benchmark standard C++ version
    print("Benchmarking C++ CPU Parity Mode...")
    t_cpp_cpu = execute_binary_repeatedly(cpp_bin, iterations, use_headers_flag=True)
    if t_cpp_cpu is None: return False

    # 3. Compile C++ binary in GPU-Hermite mode and benchmark
    print("\nRe-compiling C++ Standalone with CLOUDJ_GPU_MODE active for benchmark...")
    # Clean previous build and reconfigure with -DCLOUDJ_GPU_MODE=ON
    build_dir = os.path.dirname(os.path.dirname(fortran_bin))
    if not build_dir or build_dir == ".":
        build_dir = "build"

    try:
        subprocess.run(["cmake", "-DCLOUDJ_GPU_MODE_BENCH=ON", ".."], cwd=build_dir, check=True, stdout=subprocess.PIPE)
        subprocess.run(["make", "-j"], cwd=build_dir, check=True, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error: Re-compiling in GPU mode failed: {e}", file=sys.stderr)
        return False

    print("Benchmarking C++ GPU-Hermite Mode...")
    t_cpp_gpu = execute_binary_repeatedly(cpp_bin, iterations, use_headers_flag=True)
    if t_cpp_gpu is None: return False

    # 4. Compile C++ binary in PCR mode and benchmark
    print("\nRe-compiling C++ Standalone with CLOUDJ_USE_PCR active for benchmark...")
    try:
        subprocess.run(["cmake", "-DCLOUDJ_GPU_MODE_BENCH=OFF", "-DCLOUDJ_USE_PCR=ON", ".."], cwd=build_dir, check=True, stdout=subprocess.PIPE)
        subprocess.run(["make", "-j"], cwd=build_dir, check=True, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error: Re-compiling in PCR mode failed: {e}", file=sys.stderr)
        return False

    print("Benchmarking C++ PCR Solver Mode...")
    t_cpp_pcr = execute_binary_repeatedly(cpp_bin, iterations, use_headers_flag=True)
    if t_cpp_pcr is None: return False

    # Restore default CPU Parity mode compilation to keep repo in default state
    print("\nRestoring default CPU Parity compilation mode...")
    try:
        subprocess.run(["cmake", "-DCLOUDJ_GPU_MODE_BENCH=OFF", "-DCLOUDJ_USE_PCR=OFF", ".."], cwd=build_dir, check=True, stdout=subprocess.PIPE)
        subprocess.run(["make", "-j"], cwd=build_dir, check=True, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
         print(f"Warning: Failed to restore default compilation: {e}", file=sys.stderr)

    # 5. Compute and Print Metrics
    rate_fortran = iterations / t_fortran
    rate_cpu = iterations / t_cpp_cpu
    rate_gpu = iterations / t_cpp_gpu
    rate_pcr = iterations / t_cpp_pcr

    speedup_cpu = t_fortran / t_cpp_cpu
    speedup_gpu = t_fortran / t_cpp_gpu
    speedup_pcr = t_fortran / t_cpp_pcr
    speedup_gpu_vs_cpu = t_cpp_cpu / t_cpp_gpu

    print("\n=================================================================")
    print("                  BENCHMARKING PROFILE RESULTS                   ")
    print("=================================================================")
    print(f"Total Iterations: {iterations} vertical column scans")
    print("-----------------------------------------------------------------")
    print(f" | Compiler/Mode        | Elapsed Time | Throughput       | Speedup vs Fortran |")
    print("-----------------------------------------------------------------")
    print(f" | Fortran (gfortran)   | {t_fortran:11.4f}s | {rate_fortran:12.1f} columns/s | [Reference]        |")
    print(f" | C++ (CPU Parity)     | {t_cpp_cpu:11.4f}s | {rate_cpu:12.1f} columns/s | {speedup_cpu:17.2f}x |")
    print(f" | C++ (GPU-Hermite)    | {t_cpp_gpu:11.4f}s | {rate_gpu:12.1f} columns/s | {speedup_gpu:17.2f}x |")
    print(f" | C++ (PCR Solver)     | {t_cpp_pcr:11.4f}s | {rate_pcr:12.1f} columns/s | {speedup_pcr:17.2f}x |")
    print("=================================================================")
    print(f"Optimization Speedup (GPU-Hermite vs C++ CPU): {speedup_gpu_vs_cpu:.2f}x branchless performance increase\n")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=500)
    args = parser.parse_args()

    success = run_benchmarks(args.count)
    sys.exit(0 if success else 1)
