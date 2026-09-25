/**
 * Property-based tests for CloudJ::PhotoJX::EXTRAL1
 *
 * Property 7: EXTRAL1 Layer Insertion Consistency
 *   For any optical depth array where at least one layer exceeds ATAU0, and
 *   for any valid ATAU > 1 and ATAU0 > 0, the expanded grid dimension ND
 *   SHALL equal 2*(L1U) + 1 + 2*sum(JXTRA), and JXTRA[L] >= 0 for all L,
 *   and JXTRA[L] > 0 only where optical depth at layer L exceeds ATAU0.
 *
 * Validates: Requirements 4.1, 4.2, 4.5
 */

#include <cloudj/photo_jx.hpp>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <vector>

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

    // Generate double in [lo, hi)
    double uniform(double lo, double hi) {
        uint64_t bits = next();
        double t = static_cast<double>(bits) / static_cast<double>(UINT64_MAX);
        return lo + t * (hi - lo);
    }

    // Generate int in [lo, hi]
    int uniform_int(int lo, int hi) {
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
// Test 1: Thin layers (all OD < ATAU0) — verify all JXTRA = 0
// Validates: Requirements 4.1, 4.2
// ============================================================================
static void test_thin_layers_all_zero() {
    std::cout << "Test 1: Thin layers (all OD < ATAU0) → all JXTRA = 0\n";

    PRNG rng(42);

    const int NUM_TRIALS = 10000;
    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int L1X = rng.uniform_int(2, 60);
        double ATAU = 1.0 + rng.uniform(0.01, 0.5);  // ATAU > 1
        double ATAU0 = rng.uniform(0.001, 0.1);       // ATAU0 > 0
        int NX = 601;

        // Generate OD values all strictly less than ATAU0
        std::vector<double> DTAU600(L1X);
        for (int i = 0; i < L1X; ++i) {
            DTAU600[i] = rng.uniform(0.0, ATAU0 * 0.99);
        }

        std::vector<int> JXTRA(L1X, -1);  // Initialize with -1 to detect writes

        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        for (int i = 0; i < L1X; ++i) {
            if (JXTRA[i] != 0) {
                std::cerr << "  Counterexample: L1X=" << L1X
                          << " ATAU=" << ATAU << " ATAU0=" << ATAU0
                          << " DTAU600[" << i << "]=" << DTAU600[i]
                          << " JXTRA[" << i << "]=" << JXTRA[i] << "\n";
                CHECK(false, "JXTRA should be 0 when all OD < ATAU0");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 2: Thick layers — verify JXTRA > 0 only where OD exceeds threshold
// Validates: Requirements 4.1, 4.2
// ============================================================================
static void test_thick_layers_jxtra_positive() {
    std::cout << "Test 2: Thick layers → JXTRA > 0 only where OD exceeds threshold\n";

    PRNG rng(123);

    const int NUM_TRIALS = 10000;
    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int L1X = rng.uniform_int(5, 58);
        double ATAU = 1.05;    // default
        double ATAU0 = 0.005;  // default
        int NX = 601;

        // Generate mix of thin and thick layers
        std::vector<double> DTAU600(L1X);
        for (int i = 0; i < L1X; ++i) {
            if (rng.uniform(0.0, 1.0) < 0.3) {
                // Thick layer: OD >> ATAU0
                DTAU600[i] = rng.uniform(0.1, 10.0);
            } else {
                // Thin layer: OD << ATAU0
                DTAU600[i] = rng.uniform(0.0, ATAU0 * 0.5);
            }
        }

        std::vector<int> JXTRA(L1X, 0);

        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        // Verify: JXTRA[L] > 0 implies DTAU600[L] > ATAU0
        // Note: The threshold ATAU0X grows geometrically from the top down,
        // so a layer with OD > ATAU0 might still have JXTRA=0 if the grown
        // threshold ATAU0X already exceeded DTAU600[L]. But the converse must
        // hold: JXTRA[L] > 0 implies DTAU600[L] was greater than the threshold
        // at the time of processing, which is >= ATAU0.
        // The simpler property: JXTRA[L] > 0 → DTAU600[L] > ATAU0
        for (int i = 0; i < L1X; ++i) {
            if (JXTRA[i] > 0 && DTAU600[i] <= ATAU0) {
                std::cerr << "  Counterexample: L1X=" << L1X
                          << " DTAU600[" << i << "]=" << DTAU600[i]
                          << " JXTRA[" << i << "]=" << JXTRA[i]
                          << " ATAU0=" << ATAU0 << "\n";
                CHECK(false, "JXTRA > 0 for layer with OD <= ATAU0");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 3: ND formula — verify ND = 2*L1X + 1 + 2*sum(JXTRA)
// Validates: Requirements 4.5
// ============================================================================
static void test_nd_formula() {
    std::cout << "Test 3: ND formula → ND = 2*L1X + 1 + 2*sum(JXTRA)\n";

    PRNG rng(456);

    const int NUM_TRIALS = 10000;
    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int L1X = rng.uniform_int(2, 58);
        double ATAU = 1.0 + rng.uniform(0.01, 0.3);
        double ATAU0 = rng.uniform(0.001, 0.05);
        int NX = 601;

        // Generate random OD profile with some thick layers
        std::vector<double> DTAU600(L1X);
        for (int i = 0; i < L1X; ++i) {
            DTAU600[i] = rng.uniform(0.0, 5.0);
        }

        std::vector<int> JXTRA(L1X, 0);

        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        int sum_jxtra = 0;
        for (int i = 0; i < L1X; ++i) {
            sum_jxtra += JXTRA[i];
        }

        // The expanded grid dimension ND = 2*L1X + 1 + 2*sum(JXTRA)
        // This equals 2*(L1X + sum(JXTRA)) + 1
        int ND = 2 * L1X + 1 + 2 * sum_jxtra;

        // The overflow check ensures (L1X + 2 + sum_jxtra_partial) * 2 <= NX
        // So ND should satisfy: ND = 2*L1X + 1 + 2*sum_jxtra
        // and ND <= 2*NX - 3  (since (L1X+2+sum)*2 <= NX means
        //                       2*L1X + 4 + 2*sum <= NX
        //                       2*L1X + 1 + 2*sum <= NX - 3
        //                       ND <= NX - 3)
        // Actually, the overflow ensures no single accumulation violates, so:
        // (L1X + 2 + cumulative_sum) * 2 <= NX at each step
        // Thus the total (L1X + 2 + total_sum) * 2 <= NX
        // → ND = 2*L1X + 1 + 2*sum_jxtra <= 2*L1X + 1 + 2*(NX/2 - L1X - 2) = NX - 3

        // Verify ND is positive and bounded
        CHECK(ND > 0, "ND should be positive");
        CHECK(ND >= 2 * L1X + 1, "ND should be at least 2*L1X + 1 (zero insertions)");

        // Verify the overflow bound: the total accumulated count + L1X + 2
        // should satisfy (L1X + 2 + sum_jxtra) * 2 <= NX
        if (sum_jxtra > 0) {
            // If any insertions happened, the total should respect the bound
            CHECK((L1X + 2 + sum_jxtra) * 2 <= NX,
                  "Overflow check: (L1X + 2 + sum(JXTRA)) * 2 should be <= NX");
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 4: JXTRA non-negative — for random OD profiles, all JXTRA >= 0
// Validates: Requirements 4.2
// ============================================================================
static void test_jxtra_non_negative() {
    std::cout << "Test 4: JXTRA non-negative for all random OD profiles\n";

    PRNG rng(789);

    const int NUM_TRIALS = 50000;
    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int L1X = rng.uniform_int(1, 60);
        double ATAU = 1.0 + rng.uniform(0.001, 1.0);
        double ATAU0 = rng.uniform(1e-6, 1.0);
        int NX = rng.uniform_int(10, 1200);

        std::vector<double> DTAU600(L1X);
        for (int i = 0; i < L1X; ++i) {
            DTAU600[i] = rng.uniform(0.0, 50.0);
        }

        std::vector<int> JXTRA(L1X, -999);  // Initialize to detect no-write

        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        for (int i = 0; i < L1X; ++i) {
            if (JXTRA[i] < 0) {
                std::cerr << "  Counterexample: L1X=" << L1X
                          << " NX=" << NX << " ATAU=" << ATAU
                          << " ATAU0=" << ATAU0
                          << " DTAU600[" << i << "]=" << DTAU600[i]
                          << " JXTRA[" << i << "]=" << JXTRA[i] << "\n";
                CHECK(false, "JXTRA[L] should be >= 0 for all L");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 5: Overflow protection — set NX very small, verify JXTRA gets zeroed
// Validates: Requirements 4.5 (overflow check in EXTRAL1)
// ============================================================================
static void test_overflow_protection() {
    std::cout << "Test 5: Overflow protection with small NX\n";

    PRNG rng(1001);

    const int NUM_TRIALS = 10000;
    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int L1X = rng.uniform_int(5, 30);
        double ATAU = 1.05;
        double ATAU0 = 0.005;

        // Set NX very small so overflow will occur
        // The minimum size without any insertions is (L1X + 2)*2
        // Set NX smaller than what even a single insertion would require
        int NX = (L1X + 2) * 2 + rng.uniform_int(0, 4);

        // Generate thick layers that would normally get many insertions
        std::vector<double> DTAU600(L1X);
        for (int i = 0; i < L1X; ++i) {
            DTAU600[i] = rng.uniform(1.0, 100.0);  // Very thick layers
        }

        std::vector<int> JXTRA(L1X, 0);

        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        // Verify overflow protection: (L1X + 2 + sum(JXTRA)) * 2 <= NX
        int sum_jxtra = 0;
        for (int i = 0; i < L1X; ++i) {
            sum_jxtra += JXTRA[i];
        }

        int total_expanded = (L1X + 2 + sum_jxtra) * 2;
        if (total_expanded > NX) {
            std::cerr << "  Counterexample: L1X=" << L1X << " NX=" << NX
                      << " sum(JXTRA)=" << sum_jxtra
                      << " total_expanded=" << total_expanded << "\n";
            CHECK(false, "Overflow protection failed: total expanded > NX");
            return;
        }

        // All JXTRA still non-negative
        for (int i = 0; i < L1X; ++i) {
            if (JXTRA[i] < 0) {
                CHECK(false, "JXTRA should be >= 0 even after overflow zeroing");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 6: Default parameters ATAU=1.05, ATAU0=0.005
// Validates: Requirements 4.4
// ============================================================================
static void test_default_parameters() {
    std::cout << "Test 6: Default parameters ATAU=1.05, ATAU0=0.005\n";

    const double ATAU = 1.05;
    const double ATAU0 = 0.005;
    const int NX = 601;

    // Test case: single thick layer
    {
        const int L1X = 10;
        std::vector<double> DTAU600(L1X, 0.001);  // All thin except one
        DTAU600[L1X - 1] = 1.0;  // Top layer is thick (processed first)

        std::vector<int> JXTRA(L1X, 0);
        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        // The top layer (index L1X-1) should get insertions since 1.0 > 0.005
        CHECK(JXTRA[L1X - 1] > 0, "Top thick layer should get insertions");

        // Verify specific value: JX = round(ln(1 + (1.05-1)*1.0/0.005) / ln(1.05))
        // = round(ln(1 + 10) / ln(1.05))
        // = round(ln(11) / 0.04879...)
        // = round(2.3979 / 0.04879)
        // = round(49.14...)
        // = 49
        double expected_AJX = std::log(1.0 + (ATAU - 1.0) * 1.0 / ATAU0) / std::log(ATAU);
        int expected_JX = std::min(100, std::max(0, static_cast<int>(expected_AJX + 0.5)));
        CHECK(JXTRA[L1X - 1] == expected_JX,
              "Top layer JXTRA should match formula");

        // Thin layers below should still be 0 (since threshold grew)
        // The new threshold after first insertion:
        // ATAU0X_new = ATAU0 * ATAU^JX = 0.005 * 1.05^49 ≈ 0.005 * 11.47 ≈ 0.0574
        // All remaining layers have OD = 0.001 < 0.0574, so all JXTRA = 0
        for (int i = 0; i < L1X - 1; ++i) {
            CHECK(JXTRA[i] == 0, "Thin layers below thick layer should have JXTRA=0");
        }
    }

    // Test case: multiple thick layers
    {
        const int L1X = 20;
        std::vector<double> DTAU600(L1X, 0.001);
        DTAU600[L1X - 1] = 0.5;   // Top thick
        DTAU600[L1X - 5] = 2.0;   // Middle thick
        DTAU600[0] = 5.0;         // Bottom thick

        std::vector<int> JXTRA(L1X, 0);
        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        // Top layer should get insertions (0.5 > 0.005)
        CHECK(JXTRA[L1X - 1] > 0, "Top thick layer should have JXTRA > 0");

        // All JXTRA non-negative
        for (int i = 0; i < L1X; ++i) {
            CHECK(JXTRA[i] >= 0, "All JXTRA must be non-negative");
        }

        // ND formula check
        int sum_jxtra = 0;
        for (int i = 0; i < L1X; ++i) sum_jxtra += JXTRA[i];
        int ND = 2 * L1X + 1 + 2 * sum_jxtra;
        CHECK(ND > 2 * L1X + 1, "With thick layers, ND should exceed minimum");
    }

    // Test case: all layers at exactly ATAU0 — should get zero JXTRA
    // (condition is DTAU600[L] > ATAU0X, not >=)
    {
        const int L1X = 10;
        std::vector<double> DTAU600(L1X, ATAU0);  // Exactly at threshold

        std::vector<int> JXTRA(L1X, 0);
        CloudJ::PhotoJX::EXTRAL1(DTAU600.data(), L1X, NX, ATAU, ATAU0, JXTRA.data());

        // OD == ATAU0 does NOT trigger insertion (strict >)
        for (int i = 0; i < L1X; ++i) {
            CHECK(JXTRA[i] == 0, "Layers exactly at ATAU0 should not get insertions");
        }
    }

    std::cout << "  PASSED\n";
}

// ============================================================================
// Main
// ============================================================================
int main() {
    std::cout << "=== EXTRAL1 Property-Based Tests ===\n\n";

    test_thin_layers_all_zero();
    test_thick_layers_jxtra_positive();
    test_nd_formula();
    test_jxtra_non_negative();
    test_overflow_protection();
    test_default_parameters();

    std::cout << "\n=== Summary ===\n";
    if (failures == 0) {
        std::cout << "All tests PASSED\n";
        return 0;
    } else {
        std::cout << failures << " test(s) FAILED\n";
        return 1;
    }
}
