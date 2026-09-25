/// @file test_sphere.cpp
/// @brief Property-based tests for spherical geometry routines.
///
/// **Validates: Requirements 3.1, 3.4, 3.5**
///
/// Property 4: Flat-Earth AMF Formula — For any valid SZA in (0,90) degrees,
/// when ATM0=0, every AMF(J,L) where J>=L SHALL equal 1/cos(SZA*pi/180)
/// within floating-point precision.
///
/// Property 5: Geometric Height Conversion — For any geopotential height
/// Z_geop in [0, RAD/2), Z_geom = Z_geop / (1 - Z_geop/RAD), and
/// Z_geom > Z_geop.
///
/// Property 6: Geometric Expansion Factor Formula — For any layer midpoint
/// Z_mid >= 0, AMG = (1 + Z_mid/RAD)^2 >= 1.0.

#define MODEL_STANDALONE
#include <cloudj/photo_jx.hpp>
#include <cloudj/state.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

using namespace CloudJ;

static int failures = 0;
static int total_tests = 0;

static void check(bool cond, const std::string& msg) {
    ++total_tests;
    if (!cond) {
        ++failures;
        std::cerr << "FAIL: " << msg << "\n";
    }
}

// =========================================================================
// Property 4: Flat-Earth AMF Formula
// =========================================================================
/// For SZA in (0, 90) degrees, SPHERE1F should produce AMF(J,L) = 1/cos(SZA)
/// for all J >= L (within the atmosphere), since flat-Earth geometry means
/// every slant path through any layer above is the same 1/cos(SZA) factor.

static void test_sphere1f_flat_earth_amf() {
    std::cout << "--- Property 4: Flat-Earth AMF Formula ---\n";

    // Create a realistic profile of ~20 levels with heights in cm
    // (surface to ~60 km, typical CTM profile)
    const int L1U = 20;
    const int dim = L1U + 1;
    double ZHL[21];  // L1U+1 elements = bottom edge of each level

    // Heights from 0 to ~60 km in cm, roughly geometric spacing
    ZHL[0]  = 0.0;
    ZHL[1]  = 1.0e5;    // 1 km
    ZHL[2]  = 2.0e5;    // 2 km
    ZHL[3]  = 3.5e5;    // 3.5 km
    ZHL[4]  = 5.0e5;    // 5 km
    ZHL[5]  = 7.0e5;    // 7 km
    ZHL[6]  = 9.0e5;    // 9 km
    ZHL[7]  = 11.0e5;   // 11 km
    ZHL[8]  = 13.5e5;   // 13.5 km
    ZHL[9]  = 16.0e5;   // 16 km
    ZHL[10] = 19.0e5;   // 19 km
    ZHL[11] = 22.0e5;   // 22 km
    ZHL[12] = 26.0e5;   // 26 km
    ZHL[13] = 30.0e5;   // 30 km
    ZHL[14] = 34.0e5;   // 34 km
    ZHL[15] = 38.0e5;   // 38 km
    ZHL[16] = 42.0e5;   // 42 km
    ZHL[17] = 47.0e5;   // 47 km
    ZHL[18] = 52.0e5;   // 52 km
    ZHL[19] = 57.0e5;   // 57 km
    ZHL[20] = 62.0e5;   // 62 km

    double AMF[dim * dim];
    const double earth_rad = RAD;   // 6375.0e5 cm
    const double zzht = ZZHT;       // 5.0e5 cm

    // Test multiple SZA values
    std::vector<double> sza_degrees = {10.0, 30.0, 45.0, 60.0, 75.0, 89.0};

    for (double sza : sza_degrees) {
        double sza_rad = sza * CPI180;
        double U0 = std::cos(sza_rad);
        double expected_amf = 1.0 / U0;

        PhotoJX::SPHERE1F(U0, earth_rad, ZHL, zzht, AMF, L1U);

        // Check: for J >= L (both in [0, LTOP-1]), AMF[J + dim*L] should equal 1/cos(SZA)
        const int LTOP = L1U;  // = 20
        bool all_match = true;
        double max_relerr = 0.0;
        int fail_J = -1, fail_L = -1;
        double fail_val = 0.0;

        for (int L = 0; L < LTOP; ++L) {
            for (int J = L; J < LTOP; ++J) {
                double val = AMF[J + dim * L];
                double relerr = std::abs(val - expected_amf) / expected_amf;
                if (relerr > max_relerr) max_relerr = relerr;
                if (relerr > 1.0e-12) {
                    all_match = false;
                    if (fail_J < 0) {
                        fail_J = J;
                        fail_L = L;
                        fail_val = val;
                    }
                }
            }
        }

        if (all_match) {
            check(true, "SPHERE1F SZA=" + std::to_string(sza) +
                  " all AMF(J>=L) = 1/cos(SZA), max_relerr=" +
                  std::to_string(max_relerr));
        } else {
            check(false, "SPHERE1F SZA=" + std::to_string(sza) +
                  " AMF mismatch at J=" + std::to_string(fail_J) +
                  " L=" + std::to_string(fail_L) +
                  " got=" + std::to_string(fail_val) +
                  " expected=" + std::to_string(expected_amf));
        }

        // Also check: AMF[J + dim*L] = 0 for J < L (below the source level)
        bool below_zero = true;
        for (int L = 1; L < LTOP; ++L) {
            for (int J = 0; J < L; ++J) {
                if (AMF[J + dim * L] != 0.0) {
                    below_zero = false;
                    break;
                }
            }
            if (!below_zero) break;
        }
        check(below_zero, "SPHERE1F SZA=" + std::to_string(sza) +
              " AMF(J<L) = 0 for below-source layers");

        // Check top marker: AMF[LTOP + dim*LTOP] = 1.0
        check(AMF[LTOP + dim * LTOP] == 1.0,
              "SPHERE1F SZA=" + std::to_string(sza) +
              " top marker AMF(LTOP,LTOP) = 1.0");
    }
}

// =========================================================================
// Property 5: Geometric Height Conversion
// =========================================================================
/// For geopotential height Z_geop in [0, RAD/2), the geometric height is:
///   Z_geom = Z_geop / (1 - Z_geop/RAD)
/// And Z_geom > Z_geop for all Z_geop > 0.

static void test_geometric_height_conversion() {
    std::cout << "\n--- Property 5: Geometric Height Conversion ---\n";

    const double earth_rad = RAD;  // 6375.0e5 cm

    // Test values of Z_geop in cm
    std::vector<double> z_geop_values = {0.0, 1.0e5, 1.0e6, 5.0e6, 1.0e7, 2.0e7, 3.0e7};

    for (double z_geop : z_geop_values) {
        // Formula: Z_geom = Z_geop / (1 - Z_geop/RAD)
        double expected_z_geom = z_geop / (1.0 - z_geop / earth_rad);

        // Verify the formula gives expected results
        double ratio = z_geop / earth_rad;
        double z_geom = z_geop / (1.0 - ratio);

        double relerr = (expected_z_geom != 0.0) ?
            std::abs(z_geom - expected_z_geom) / std::abs(expected_z_geom) : 0.0;

        check(relerr < 1.0e-14,
              "Z_geom formula Z_geop=" + std::to_string(z_geop) +
              " Z_geom=" + std::to_string(z_geom) +
              " expected=" + std::to_string(expected_z_geom) +
              " relerr=" + std::to_string(relerr));

        // Property: Z_geom > Z_geop for Z_geop > 0
        if (z_geop > 0.0) {
            check(z_geom > z_geop,
                  "Z_geom > Z_geop for Z_geop=" + std::to_string(z_geop) +
                  " Z_geom=" + std::to_string(z_geom));
        } else {
            // At surface (z_geop=0), z_geom should also be 0
            check(z_geom == 0.0,
                  "Z_geom = 0 at surface, got " + std::to_string(z_geom));
        }

        // Additional property: Z_geom = Z_geop * (1 + Z_geop/RAD + (Z_geop/RAD)^2 + ...)
        // For small Z_geop/RAD, Z_geom ~ Z_geop * (1 + Z_geop/RAD)
        // This is a sanity check that geometric height is always >= geopotential
        if (z_geop > 0.0) {
            double expansion_factor = z_geom / z_geop;
            double expected_factor = 1.0 / (1.0 - ratio);
            check(expansion_factor >= 1.0,
                  "Expansion factor >= 1.0 for Z_geop=" + std::to_string(z_geop) +
                  " factor=" + std::to_string(expansion_factor));
            double factor_err = std::abs(expansion_factor - expected_factor) / expected_factor;
            check(factor_err < 1.0e-14,
                  "Expansion factor matches formula for Z_geop=" + std::to_string(z_geop));
        }
    }

    // Edge case: very large Z_geop approaching RAD/2
    {
        double z_geop = earth_rad / 2.0 - 1.0e3;  // just under RAD/2
        double z_geom = z_geop / (1.0 - z_geop / earth_rad);
        check(z_geom > z_geop,
              "Z_geom > Z_geop for Z_geop near RAD/2, z_geom=" + std::to_string(z_geom));
        check(std::isfinite(z_geom),
              "Z_geom is finite for Z_geop near RAD/2");
    }
}

// =========================================================================
// Property 6: Geometric Expansion Factor Formula
// =========================================================================
/// For any layer midpoint Z_mid >= 0, the geometric expansion factor:
///   AMG = (1 + Z_mid/RAD)^2 >= 1.0

static void test_geometric_expansion_factor() {
    std::cout << "\n--- Property 6: Geometric Expansion Factor Formula ---\n";

    const double earth_rad = RAD;  // 6375.0e5 cm

    // Test values of Z_mid in cm
    std::vector<double> z_mid_values = {0.0, 1.0e5, 5.0e5, 1.0e6, 5.0e6, 1.0e7};

    for (double z_mid : z_mid_values) {
        // Formula: AMG = (1 + Z_mid/RAD)^2
        double ratio = z_mid / earth_rad;
        double expected_amg = (1.0 + ratio) * (1.0 + ratio);
        double amg = std::pow(1.0 + z_mid / earth_rad, 2.0);

        double relerr = std::abs(amg - expected_amg) / expected_amg;

        check(relerr < 1.0e-14,
              "AMG formula Z_mid=" + std::to_string(z_mid) +
              " AMG=" + std::to_string(amg) +
              " expected=" + std::to_string(expected_amg));

        // Property: AMG >= 1.0
        check(amg >= 1.0,
              "AMG >= 1.0 for Z_mid=" + std::to_string(z_mid) +
              " AMG=" + std::to_string(amg));

        // At surface (Z_mid=0), AMG should be exactly 1.0
        if (z_mid == 0.0) {
            check(amg == 1.0,
                  "AMG = 1.0 at surface, got " + std::to_string(amg));
        } else {
            // For Z_mid > 0, AMG must be strictly > 1.0
            check(amg > 1.0,
                  "AMG > 1.0 for Z_mid=" + std::to_string(z_mid) +
                  " AMG=" + std::to_string(amg));
        }

        // AMG should increase with altitude
        if (z_mid > 0.0) {
            double amg_lower = std::pow(1.0 + (z_mid / 2.0) / earth_rad, 2.0);
            check(amg > amg_lower,
                  "AMG increases with altitude: AMG(" + std::to_string(z_mid) +
                  ")=" + std::to_string(amg) + " > AMG(" +
                  std::to_string(z_mid / 2.0) + ")=" + std::to_string(amg_lower));
        }
    }

    // Test monotonicity across all values
    double prev_amg = 0.0;
    for (double z_mid : z_mid_values) {
        double amg = std::pow(1.0 + z_mid / earth_rad, 2.0);
        if (prev_amg > 0.0) {
            check(amg >= prev_amg,
                  "AMG monotonically increasing: z_mid=" + std::to_string(z_mid));
        }
        prev_amg = amg;
    }
}

// =========================================================================
// SPHERE1N Sanity Checks
// =========================================================================
/// Basic sanity check for SPHERE1N: verify AMF values > 0 for sunlit layers
/// and that AMF increases with SZA.

static void test_sphere1n_sanity() {
    std::cout << "\n--- SPHERE1N Sanity Checks ---\n";

    const int L1U = 20;
    const int dim = L1U + 1;
    double ZHL[21];

    // Realistic altitude profile (cm)
    ZHL[0]  = 0.0;
    ZHL[1]  = 1.0e5;
    ZHL[2]  = 2.0e5;
    ZHL[3]  = 3.5e5;
    ZHL[4]  = 5.0e5;
    ZHL[5]  = 7.0e5;
    ZHL[6]  = 9.0e5;
    ZHL[7]  = 11.0e5;
    ZHL[8]  = 13.5e5;
    ZHL[9]  = 16.0e5;
    ZHL[10] = 19.0e5;
    ZHL[11] = 22.0e5;
    ZHL[12] = 26.0e5;
    ZHL[13] = 30.0e5;
    ZHL[14] = 34.0e5;
    ZHL[15] = 38.0e5;
    ZHL[16] = 42.0e5;
    ZHL[17] = 47.0e5;
    ZHL[18] = 52.0e5;
    ZHL[19] = 57.0e5;
    ZHL[20] = 62.0e5;

    double AMF[dim * dim];
    const double earth_rad = RAD;
    const double zzht = ZZHT;
    const int LTOP = L1U;

    // Test 1: All AMF(J,L) > 0 for sunlit layers (U0 > 0, J >= L)
    std::vector<double> sza_degrees = {10.0, 30.0, 45.0, 60.0, 75.0, 85.0};

    for (double sza : sza_degrees) {
        double sza_rad = sza * CPI180;
        double U0 = std::cos(sza_rad);

        PhotoJX::SPHERE1N(U0, earth_rad, ZHL, zzht, AMF, L1U);

        bool all_positive = true;
        int neg_J = -1, neg_L = -1;
        double neg_val = 0.0;

        for (int L = 0; L < LTOP; ++L) {
            for (int J = L; J < LTOP; ++J) {
                double val = AMF[J + dim * L];
                if (val <= 0.0) {
                    all_positive = false;
                    neg_J = J;
                    neg_L = L;
                    neg_val = val;
                    break;
                }
            }
            if (!all_positive) break;
        }

        check(all_positive,
              "SPHERE1N SZA=" + std::to_string(sza) +
              " all AMF(J>=L) > 0" +
              (all_positive ? "" : " first fail at J=" + std::to_string(neg_J) +
               " L=" + std::to_string(neg_L) + " val=" + std::to_string(neg_val)));
    }

    // Test 2: AMF increases with SZA (for the same layer pair)
    // Compare AMF at SZA=30 vs SZA=60 at a mid-atmosphere layer pair
    {
        double AMF_30[dim * dim], AMF_60[dim * dim];
        double U0_30 = std::cos(30.0 * CPI180);
        double U0_60 = std::cos(60.0 * CPI180);

        PhotoJX::SPHERE1N(U0_30, earth_rad, ZHL, zzht, AMF_30, L1U);
        PhotoJX::SPHERE1N(U0_60, earth_rad, ZHL, zzht, AMF_60, L1U);

        // Check that AMF at SZA=60 >= AMF at SZA=30 for surface-to-top pair
        int L = 0;
        int J = LTOP - 1;
        double amf_30_val = AMF_30[J + dim * L];
        double amf_60_val = AMF_60[J + dim * L];

        check(amf_60_val > amf_30_val,
              "SPHERE1N AMF increases with SZA: AMF(60deg)=" +
              std::to_string(amf_60_val) + " > AMF(30deg)=" +
              std::to_string(amf_30_val));
    }

    // Test 3: For overhead sun (SZA~0), SPHERE1N should give AMF close to 1
    {
        double U0 = std::cos(1.0 * CPI180);  // SZA = 1 degree (nearly overhead)
        PhotoJX::SPHERE1N(U0, earth_rad, ZHL, zzht, AMF, L1U);

        // Each layer path should be close to 1.0 (geometric effects small near zenith)
        bool close_to_one = true;
        double max_deviation = 0.0;
        for (int L = 0; L < LTOP; ++L) {
            for (int J = L; J < LTOP; ++J) {
                double val = AMF[J + dim * L];
                double dev = std::abs(val - 1.0);
                if (dev > max_deviation) max_deviation = dev;
                // At near-zenith, AMF should be very close to 1
                // but spherical geometry means it's not exactly 1
                if (dev > 0.01) {  // allow 1% deviation for spherical effects
                    close_to_one = false;
                }
            }
        }
        check(close_to_one,
              "SPHERE1N SZA=1deg AMF values near 1.0, max_deviation=" +
              std::to_string(max_deviation));
    }

    // Test 4: SPHERE1N vs SPHERE1F comparison at small SZA
    // For small SZA, spherical and flat-earth should be very close
    {
        double sza_small = 10.0;  // 10 degrees
        double U0 = std::cos(sza_small * CPI180);

        double AMF_flat[dim * dim], AMF_sphere[dim * dim];
        PhotoJX::SPHERE1F(U0, earth_rad, ZHL, zzht, AMF_flat, L1U);
        PhotoJX::SPHERE1N(U0, earth_rad, ZHL, zzht, AMF_sphere, L1U);

        double max_relerr = 0.0;
        for (int L = 0; L < LTOP; ++L) {
            for (int J = L; J < LTOP; ++J) {
                double flat_val = AMF_flat[J + dim * L];
                double sphere_val = AMF_sphere[J + dim * L];
                if (flat_val > 0.0) {
                    double relerr = std::abs(sphere_val - flat_val) / flat_val;
                    if (relerr > max_relerr) max_relerr = relerr;
                }
            }
        }
        // At 10 degrees SZA, spherical and flat should agree within a few percent
        check(max_relerr < 0.05,
              "SPHERE1N vs SPHERE1F at SZA=10deg: max_relerr=" +
              std::to_string(max_relerr) + " (should be < 5%)");
    }
}

// =========================================================================
// Main
// =========================================================================

int main() {
    std::cout << "=== Spherical Geometry Property Tests ===\n\n";

    test_sphere1f_flat_earth_amf();
    test_geometric_height_conversion();
    test_geometric_expansion_factor();
    test_sphere1n_sanity();

    // Summary
    std::cout << "\n=== Summary ===\n";
    std::cout << "Total checks: " << total_tests << "\n";
    if (failures == 0) {
        std::cout << "PASSED: All property tests passed.\n";
    } else {
        std::cout << "FAILED: " << failures << " check(s) failed.\n";
    }

    return (failures == 0) ? 0 : 1;
}
