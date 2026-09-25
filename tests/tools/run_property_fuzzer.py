# tests/tools/run_property_fuzzer.py
import sys
import os
import subprocess
import random
import argparse

def generate_random_math_input():
    # Temperature bounded between 150K and 450K
    temp = random.uniform(150.0, 450.0)
    # Target bounds
    t1 = random.uniform(180.0, 220.0)
    t2 = random.uniform(280.0, 320.0)
    t3 = random.uniform(380.0, 420.0)
    # Cross-section variables
    x1 = random.uniform(1e-22, 1e-18)
    x2 = random.uniform(1e-22, 1e-18)
    x3 = random.uniform(1e-22, 1e-18)
    # Number of interpolation points (1, 2, or 3)
    lqq = random.choice([1, 2, 3])
    
    return f"{temp:.9f} {t1:.9f} {x1:.9E} {t2:.9f} {x2:.9E} {t3:.9f} {x3:.9E} {lqq}\n"

def execute_fuzzer(count):
    print(f"Running randomized physical property fuzzer across {count} iterations...")
    
    fortran_bin = "./test_fortran_math"
    cpp_bin = "./test_cpp_math"
    
    # Locate paths dynamically
    if not os.path.exists(fortran_bin) or not os.path.exists(cpp_bin):
        for path in ["build/bin", "bin", "."]:
            f_path = os.path.join(path, "test_fortran_math")
            c_path = os.path.join(path, "test_cpp_math")
            if os.path.exists(f_path) and os.path.exists(c_path):
                fortran_bin = f_path
                cpp_bin = c_path
                break
        else:
            print("Error: Standalone math binaries not found in current directory.", file=sys.stderr)
            return False

    print(f"Running math comparison binaries located at:\n  Fortran: {fortran_bin}\n  C++:     {cpp_bin}\n")

    # Spawn both subprocesses
    proc_f = subprocess.Popen([fortran_bin], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    proc_c = subprocess.Popen([cpp_bin], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)

    disparities = 0
    tolerance = 1e-12

    for i in range(1, count + 1):
        input_data = generate_random_math_input()
        
        # Write to processes
        proc_f.stdin.write(input_data)
        proc_f.stdin.flush()
        proc_c.stdin.write(input_data)
        proc_c.stdin.flush()

        # Read results
        res_f_str = proc_f.stdout.readline().strip()
        res_c_str = proc_c.stdout.readline().strip()

        if not res_f_str or not res_c_str:
            print("Error: Subprocess returned empty stream.", file=sys.stderr)
            break

        val_f = float(res_f_str)
        val_c = float(res_c_str)

        # Compute delta differences
        diff = abs(val_f - val_c)
        rel_diff = diff / max(1e-30, abs(val_f))

        # Check precision epsilon
        if diff > tolerance and rel_diff > tolerance:
            print(f"Disparity found at iteration {i}:", file=sys.stderr)
            print(f"  Inputs:  {input_data.strip()}", file=sys.stderr)
            print(f"  Fortran: {val_f:.15E}", file=sys.stderr)
            print(f"  C++:     {val_c:.15E}", file=sys.stderr)
            print(f"  Diff:    {diff:.15E}", file=sys.stderr)
            disparities += 1
            break

        if i % 1000 == 0:
            print(f"Successfully processed {i} / {count} randomized physical calculations...")

    # Close processes
    proc_f.stdin.close()
    proc_c.stdin.close()
    proc_f.wait()
    proc_c.wait()

    if disparities > 0:
        print("Fuzz test verification failed!", file=sys.stderr)
        return False

    print("=================================================================")
    print(f"SUCCESS: Physical property fuzzer passed all {count} iterations!")
    print("=================================================================")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=10000)
    args = parser.parse_args()
    
    success = execute_fuzzer(args.count)
    sys.exit(0 if success else 1)
