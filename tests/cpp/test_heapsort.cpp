/**
 * Property-based tests for CloudJ::CloudOverlap::HEAPSORT_A
 *
 * Property 12: HEAPSORT Produces Sorted Output
 *   For any array of N positive doubles, after calling HEAPSORT_A, the output
 *   array SHALL be sorted in non-decreasing order and contain exactly the same
 *   elements as the input (permutation invariant).
 *
 * Validates: Requirements 7.10
 */

#include <cloudj/cloud_jx.hpp>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstring>

// Simple xorshift64 PRNG for property-based testing
struct PRNG {
    uint64_t state;

    explicit PRNG(uint64_t seed) : state(seed ? seed : 1) {}

    uint64_t next() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }

    // Generate a random double in (0, max_val]
    double random_positive_double(double max_val = 1e6) {
        uint64_t bits = next();
        double t = static_cast<double>(bits) / static_cast<double>(UINT64_MAX);
        return t * max_val + 1e-15; // ensure strictly positive
    }

    // Generate a random int in [lo, hi]
    int random_int(int lo, int hi) {
        return lo + static_cast<int>(next() % static_cast<uint64_t>(hi - lo + 1));
    }
};

static int failures = 0;

static void check(bool condition, const char* msg, int line) {
    if (!condition) {
        std::cerr << "FAIL (line " << line << "): " << msg << "\n";
        ++failures;
    }
}

#define CHECK(cond, msg) check((cond), (msg), __LINE__)

// ============================================================================
// Helper: Check that output array is sorted in non-decreasing order
// ============================================================================
static bool is_sorted_nondecreasing(const double* AX, int N) {
    for (int i = 1; i < N; ++i) {
        if (AX[i] < AX[i - 1]) return false;
    }
    return true;
}

// ============================================================================
// Helper: Check that sorted array is a permutation of the input
// (same multiset of elements)
// ============================================================================
static bool is_permutation_of(const double* A, const double* AX, int N) {
    std::vector<double> input(A, A + N);
    std::vector<double> output(AX, AX + N);
    std::sort(input.begin(), input.end());
    std::sort(output.begin(), output.end());
    return input == output;
}

// ============================================================================
// Helper: Check that IX permutation correctly maps original to sorted positions
// i.e., A[IX[j]] == AX[j] for all j
// ============================================================================
static bool index_permutation_valid(const double* A, const double* AX, const int* IX, int N) {
    for (int j = 0; j < N; ++j) {
        if (IX[j] < 0 || IX[j] >= N) return false;
        if (A[IX[j]] != AX[j]) return false;
    }
    return true;
}

// ============================================================================
// Property 12: HEAPSORT Produces Sorted Output
// Validates: Requirements 7.10
//
// Sub-properties tested:
//   12a: Output is sorted in non-decreasing order
//   12b: Output is a permutation of the input (same elements)
//   12c: Index array IX correctly maps original positions to sorted positions
// ============================================================================

// --- Test 1: Random arrays ---
static void test_random_arrays() {
    std::cout << "  Random arrays (N=2..100, 10000 trials)...\n";

    PRNG rng(42);
    const int NUM_TRIALS = 10000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int N = rng.random_int(2, 100);
        int ND = N; // ND >= N

        std::vector<double> A(N);
        std::vector<double> AX(N);
        std::vector<int> IX(N);

        // Generate random positive doubles
        for (int i = 0; i < N; ++i) {
            A[i] = rng.random_positive_double();
        }

        CloudJ::CloudOverlap::HEAPSORT_A(N, A.data(), AX.data(), IX.data(), ND);

        // Property 12a: sorted non-decreasing
        if (!is_sorted_nondecreasing(AX.data(), N)) {
            std::cerr << "  Counterexample (trial " << trial << "): N=" << N
                      << " output not sorted\n";
            std::cerr << "  Input:  ";
            for (int i = 0; i < std::min(N, 10); ++i) std::cerr << A[i] << " ";
            std::cerr << "\n  Output: ";
            for (int i = 0; i < std::min(N, 10); ++i) std::cerr << AX[i] << " ";
            std::cerr << "\n";
            CHECK(false, "HEAPSORT_A output not sorted (random array)");
            return;
        }

        // Property 12b: permutation invariant
        if (!is_permutation_of(A.data(), AX.data(), N)) {
            std::cerr << "  Counterexample (trial " << trial << "): N=" << N
                      << " output is not a permutation of input\n";
            CHECK(false, "HEAPSORT_A output not a permutation of input (random array)");
            return;
        }

        // Property 12c: index permutation valid
        if (!index_permutation_valid(A.data(), AX.data(), IX.data(), N)) {
            std::cerr << "  Counterexample (trial " << trial << "): N=" << N
                      << " IX permutation invalid\n";
            CHECK(false, "HEAPSORT_A IX permutation invalid (random array)");
            return;
        }
    }

    std::cout << "    PASSED\n";
}

// --- Test 2: Permutation check with duplicate values ---
static void test_permutation_with_duplicates() {
    std::cout << "  Permutation check with duplicate values...\n";

    PRNG rng(99);
    const int NUM_TRIALS = 5000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int N = rng.random_int(2, 50);
        int ND = N;

        std::vector<double> A(N);
        std::vector<double> AX(N);
        std::vector<int> IX(N);

        // Generate array with many duplicate values (small value range)
        int num_distinct = rng.random_int(1, 5);
        std::vector<double> vals(num_distinct);
        for (int i = 0; i < num_distinct; ++i) {
            vals[i] = rng.random_positive_double(10.0);
        }
        for (int i = 0; i < N; ++i) {
            A[i] = vals[rng.random_int(0, num_distinct - 1)];
        }

        CloudJ::CloudOverlap::HEAPSORT_A(N, A.data(), AX.data(), IX.data(), ND);

        if (!is_sorted_nondecreasing(AX.data(), N)) {
            std::cerr << "  Counterexample (trial " << trial << "): N=" << N
                      << " output not sorted with duplicates\n";
            CHECK(false, "HEAPSORT_A output not sorted (duplicates)");
            return;
        }

        if (!is_permutation_of(A.data(), AX.data(), N)) {
            std::cerr << "  Counterexample (trial " << trial << "): N=" << N
                      << " output not a permutation (duplicates)\n";
            CHECK(false, "HEAPSORT_A output not permutation (duplicates)");
            return;
        }

        if (!index_permutation_valid(A.data(), AX.data(), IX.data(), N)) {
            std::cerr << "  Counterexample (trial " << trial << "): N=" << N
                      << " IX invalid (duplicates)\n";
            CHECK(false, "HEAPSORT_A IX invalid (duplicates)");
            return;
        }
    }

    std::cout << "    PASSED\n";
}

// --- Test 3: Edge cases ---
static void test_edge_cases() {
    std::cout << "  Edge cases (N=0, N=1, all-same, already-sorted, reverse-sorted)...\n";

    // N=0: should not crash
    {
        std::vector<double> AX(1, -1.0);
        std::vector<int> IX(1, -1);
        CloudJ::CloudOverlap::HEAPSORT_A(0, nullptr, AX.data(), IX.data(), 0);
        // Nothing to verify except no crash
        std::cout << "    N=0: OK (no crash)\n";
    }

    // N=1: single element
    {
        double A[1] = {3.14};
        double AX[1] = {0.0};
        int IX[1] = {-1};
        CloudJ::CloudOverlap::HEAPSORT_A(1, A, AX, IX, 1);
        CHECK(AX[0] == 3.14, "N=1: AX[0] should equal input");
        CHECK(IX[0] == 0, "N=1: IX[0] should be 0");
        std::cout << "    N=1: OK\n";
    }

    // All-same values
    {
        const int N = 20;
        double A[N];
        double AX[N];
        int IX[N];
        for (int i = 0; i < N; ++i) A[i] = 7.7;

        CloudJ::CloudOverlap::HEAPSORT_A(N, A, AX, IX, N);

        CHECK(is_sorted_nondecreasing(AX, N), "All-same: should be sorted");
        CHECK(is_permutation_of(A, AX, N), "All-same: should be permutation");
        CHECK(index_permutation_valid(A, AX, IX, N), "All-same: IX should be valid");
        std::cout << "    All-same (N=20): OK\n";
    }

    // Already sorted
    {
        const int N = 50;
        double A[N];
        double AX[N];
        int IX[N];
        for (int i = 0; i < N; ++i) A[i] = static_cast<double>(i + 1);

        CloudJ::CloudOverlap::HEAPSORT_A(N, A, AX, IX, N);

        CHECK(is_sorted_nondecreasing(AX, N), "Already-sorted: should remain sorted");
        CHECK(is_permutation_of(A, AX, N), "Already-sorted: should be permutation");
        CHECK(index_permutation_valid(A, AX, IX, N), "Already-sorted: IX should be valid");
        // Verify identity permutation (since input is already sorted)
        bool identity = true;
        for (int i = 0; i < N; ++i) {
            if (IX[i] != i) { identity = false; break; }
        }
        CHECK(identity, "Already-sorted: IX should be identity permutation");
        std::cout << "    Already-sorted (N=50): OK\n";
    }

    // Reverse-sorted
    {
        const int N = 50;
        double A[N];
        double AX[N];
        int IX[N];
        for (int i = 0; i < N; ++i) A[i] = static_cast<double>(N - i);

        CloudJ::CloudOverlap::HEAPSORT_A(N, A, AX, IX, N);

        CHECK(is_sorted_nondecreasing(AX, N), "Reverse-sorted: should be sorted ascending");
        CHECK(is_permutation_of(A, AX, N), "Reverse-sorted: should be permutation");
        CHECK(index_permutation_valid(A, AX, IX, N), "Reverse-sorted: IX should be valid");
        // Verify reverse permutation: IX[i] should be N-1-i
        bool reversed = true;
        for (int i = 0; i < N; ++i) {
            if (IX[i] != N - 1 - i) { reversed = false; break; }
        }
        CHECK(reversed, "Reverse-sorted: IX should map reversed positions");
        std::cout << "    Reverse-sorted (N=50): OK\n";
    }

    std::cout << "    PASSED\n";
}

// --- Test 4: Large arrays ---
static void test_large_arrays() {
    std::cout << "  Large arrays (N=100, N=1000)...\n";

    PRNG rng(777);

    for (int N : {100, 1000}) {
        int ND = N;
        std::vector<double> A(N);
        std::vector<double> AX(N);
        std::vector<int> IX(N);

        for (int i = 0; i < N; ++i) {
            A[i] = rng.random_positive_double(1e8);
        }

        CloudJ::CloudOverlap::HEAPSORT_A(N, A.data(), AX.data(), IX.data(), ND);

        CHECK(is_sorted_nondecreasing(AX.data(), N),
              (std::string("N=") + std::to_string(N) + ": output not sorted").c_str());
        CHECK(is_permutation_of(A.data(), AX.data(), N),
              (std::string("N=") + std::to_string(N) + ": output not permutation").c_str());
        CHECK(index_permutation_valid(A.data(), AX.data(), IX.data(), N),
              (std::string("N=") + std::to_string(N) + ": IX invalid").c_str());

        std::cout << "    N=" << N << ": OK\n";
    }

    std::cout << "    PASSED\n";
}

// --- Test 5: Input array is not modified (const correctness) ---
static void test_input_not_modified() {
    std::cout << "  Input array not modified (const correctness)...\n";

    PRNG rng(2023);
    const int N = 30;
    int ND = N;

    std::vector<double> A(N);
    std::vector<double> A_copy(N);
    std::vector<double> AX(N);
    std::vector<int> IX(N);

    for (int i = 0; i < N; ++i) {
        A[i] = rng.random_positive_double();
    }
    A_copy = A;

    CloudJ::CloudOverlap::HEAPSORT_A(N, A.data(), AX.data(), IX.data(), ND);

    // Verify input A is unchanged
    bool unchanged = true;
    for (int i = 0; i < N; ++i) {
        if (A[i] != A_copy[i]) { unchanged = false; break; }
    }
    CHECK(unchanged, "Input array A should not be modified by HEAPSORT_A");

    std::cout << "    PASSED\n";
}

// --- Test 6: IX contains a valid permutation of indices [0, N) ---
static void test_ix_is_valid_permutation_of_indices() {
    std::cout << "  IX is a valid permutation of indices [0, N)...\n";

    PRNG rng(314);
    const int NUM_TRIALS = 5000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int N = rng.random_int(2, 80);
        int ND = N;

        std::vector<double> A(N);
        std::vector<double> AX(N);
        std::vector<int> IX(N);

        for (int i = 0; i < N; ++i) {
            A[i] = rng.random_positive_double(100.0);
        }

        CloudJ::CloudOverlap::HEAPSORT_A(N, A.data(), AX.data(), IX.data(), ND);

        // Check that IX contains each index 0..N-1 exactly once
        std::vector<int> sorted_ix(IX.begin(), IX.begin() + N);
        std::sort(sorted_ix.begin(), sorted_ix.end());
        bool valid = true;
        for (int i = 0; i < N; ++i) {
            if (sorted_ix[i] != i) { valid = false; break; }
        }

        if (!valid) {
            std::cerr << "  Counterexample (trial " << trial << "): N=" << N
                      << " IX not a valid permutation of [0,N)\n";
            CHECK(false, "HEAPSORT_A IX not a valid permutation of indices");
            return;
        }
    }

    std::cout << "    PASSED\n";
}

int main() {
    std::cout << "=== HEAPSORT_A Property-Based Tests ===\n";
    std::cout << "Property 12: HEAPSORT Produces Sorted Output\n";
    std::cout << "Validates: Requirements 7.10\n\n";

    test_random_arrays();
    test_permutation_with_duplicates();
    test_edge_cases();
    test_large_arrays();
    test_input_not_modified();
    test_ix_is_valid_permutation_of_indices();

    std::cout << "\n=== Summary ===\n";
    if (failures == 0) {
        std::cout << "All tests PASSED\n";
        return 0;
    } else {
        std::cout << failures << " test(s) FAILED\n";
        return 1;
    }
}
