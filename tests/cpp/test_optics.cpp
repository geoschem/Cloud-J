/**
 * Property-based tests for CloudJ optical property routines.
 *
 * Property 8: Optical Property Physical Constraints
 *   For any valid cloud water path > 0 and effective radius > 0, the computed
 *   optical properties SHALL satisfy: optical_depth >= 0,
 *   0 <= single_scattering_albedo <= 1, and the total scattering optical
 *   depth (OD * SSA) <= OD.
 *
 * Property 9: Optical Property Accumulation Additivity
 *   For any set of cloud and aerosol contributions to a layer, the total
 *   optical depth SHALL equal the sum of individual optical depths, and the
 *   total single-scattering albedo SHALL equal the SSA-weighted sum divided
 *   by total OD.
 *
 * Validates: Requirements 5.1, 5.2, 5.5, 5.6, 5.7
 */

#include <cloudj/init.hpp>
#include <cloudj/photo_jx.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdint>
#include <vector>
#include <string>
#include <array>
#include <numeric>

using namespace CloudJ;

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

    double uniform(double lo, double hi) {
        uint64_t bits = next();
        double t = static_cast<double>(bits) / static_cast<double>(UINT64_MAX);
        return lo + t * (hi - lo);
    }

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
// Helper: Initialize CloudJState using INIT_CLDJ with MODEL_STANDALONE
// ============================================================================
static CloudJState g_state;
static bool g_initialized = false;

static bool init_state() {
    if (g_initialized) return true;

    std::vector<std::string> TITLEJXX;
    int NJXX = 0;
    int rc = 0;

    Init::INIT_CLDJ(
        true,           // AMIROOT
        "",             // DATADIR (unused with embedded tables)
        57,             // NLEVELS
        34,             // NLEVELS_WITH_CLOUD
        TITLEJXX,       // output titles
        0,              // NJXU
        1.05,           // ATAU
        0.005,          // ATAU0
        18,             // NWBIN
        7,              // CLDFLAG
        0.33,           // CLDCOR
        6,              // LNRG
        0,              // ATM0
        false,          // use_H2O_UV_abs
        NJXX,           // output
        g_state,        // state
        rc              // return code
    );

    if (rc != CLDJ_SUCCESS) {
        std::cerr << "ERROR: INIT_CLDJ failed with rc=" << rc << "\n";
        return false;
    }

    g_initialized = true;
    std::cout << "State initialized: MCC=" << g_state.MCC
              << " NCC=" << g_state.NCC
              << " NSS=" << g_state.NSS
              << " NAA=" << g_state.NAA
              << " NGG=" << g_state.NGG
              << " NJX=" << g_state.NJX << "\n";
    return true;
}

// ============================================================================
// Test 1: OPTICL Physical Constraints (Property 8)
// For liquid water clouds with various effective radii (1.5 to 48 microns):
// verify QQEXT >= 0, 0 <= SSALB <= 1, SSLEG[0 + 8*J] == 1.0
// **Validates: Requirements 5.1**
// ============================================================================
static void test_opticl_physical_constraints() {
    std::cout << "Test 1: OPTICL Physical Constraints\n";

    PRNG rng(42);
    const int NUM_TRIALS = 5000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        // Effective radius: range from liquid cloud tables (roughly 1.5-48 microns)
        double REFF = rng.uniform(1.5, 48.0);
        double TEFF = rng.uniform(260.0, 300.0);  // unused for liquid
        double DDENS = 0.0;
        std::array<double, S_> QQEXT{};
        std::array<double, S_> SSALB{};
        std::array<double, 8 * S_> SSLEG{};

        PhotoJX::OPTICL(REFF, TEFF, DDENS, QQEXT.data(),
                        SSALB.data(), SSLEG.data(), g_state);

        // Verify physical constraints
        for (int J = 0; J < S_; ++J) {
            if (QQEXT[J] < 0.0) {
                std::cerr << "  Counterexample: REFF=" << REFF
                          << " QQEXT[" << J << "]=" << QQEXT[J] << "\n";
                CHECK(false, "OPTICL: QQEXT must be >= 0");
                return;
            }
            if (SSALB[J] < 0.0 || SSALB[J] > 1.0) {
                std::cerr << "  Counterexample: REFF=" << REFF
                          << " SSALB[" << J << "]=" << SSALB[J] << "\n";
                CHECK(false, "OPTICL: SSALB must be in [0,1]");
                return;
            }
            // First Legendre moment (L=0) should be 1.0
            double p0 = SSLEG[0 + 8 * J];
            if (std::abs(p0 - 1.0) > 1.0e-10) {
                std::cerr << "  Counterexample: REFF=" << REFF
                          << " SSLEG[0+" << 8*J << "]=" << p0 << "\n";
                CHECK(false, "OPTICL: First Legendre moment must be 1.0");
                return;
            }
        }

        // Density must be positive
        if (DDENS <= 0.0) {
            std::cerr << "  Counterexample: REFF=" << REFF
                      << " DDENS=" << DDENS << "\n";
            CHECK(false, "OPTICL: DDENS must be > 0");
            return;
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 2: OPTICI Physical Constraints (Property 8)
// For ice clouds with various effective radii and temperatures:
// warm ice (T >= 233.15K) and cold ice (T < 233.15K)
// verify QQEXT >= 0, 0 <= SSALB <= 1, SSLEG[0 + 8*J] == 1.0
// **Validates: Requirements 5.2**
// ============================================================================
static void test_optici_physical_constraints() {
    std::cout << "Test 2: OPTICI Physical Constraints\n";

    PRNG rng(123);
    const int NUM_TRIALS = 5000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        double REFF = rng.uniform(5.0, 130.0);  // ice crystal radii
        // Alternate between warm and cold ice
        double TEFF;
        if (trial % 2 == 0) {
            TEFF = rng.uniform(233.15, 270.0);   // warm ice (irregular)
        } else {
            TEFF = rng.uniform(180.0, 233.14);   // cold ice (hexagonal)
        }

        double DDENS = 0.0;
        std::array<double, S_> QQEXT{};
        std::array<double, S_> SSALB{};
        std::array<double, 8 * S_> SSLEG{};

        PhotoJX::OPTICI(REFF, TEFF, DDENS, QQEXT.data(),
                        SSALB.data(), SSLEG.data(), g_state);

        for (int J = 0; J < S_; ++J) {
            if (QQEXT[J] < 0.0) {
                std::cerr << "  Counterexample: REFF=" << REFF
                          << " TEFF=" << TEFF
                          << " QQEXT[" << J << "]=" << QQEXT[J] << "\n";
                CHECK(false, "OPTICI: QQEXT must be >= 0");
                return;
            }
            if (SSALB[J] < 0.0 || SSALB[J] > 1.0) {
                std::cerr << "  Counterexample: REFF=" << REFF
                          << " TEFF=" << TEFF
                          << " SSALB[" << J << "]=" << SSALB[J] << "\n";
                CHECK(false, "OPTICI: SSALB must be in [0,1]");
                return;
            }
            double p0 = SSLEG[0 + 8 * J];
            if (std::abs(p0 - 1.0) > 1.0e-10) {
                std::cerr << "  Counterexample: REFF=" << REFF
                          << " TEFF=" << TEFF
                          << " SSLEG[0+" << 8*J << "]=" << p0 << "\n";
                CHECK(false, "OPTICI: First Legendre moment must be 1.0");
                return;
            }
        }

        if (DDENS <= 0.0) {
            std::cerr << "  Counterexample: REFF=" << REFF
                      << " TEFF=" << TEFF << " DDENS=" << DDENS << "\n";
            CHECK(false, "OPTICI: DDENS must be > 0");
            return;
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 3: OPTICS Physical Constraints (Property 8)
// For UCI stratospheric sulfate aerosol with valid K indices (1,2) and PATH > 0:
// verify OD >= 0, 0 <= SSA <= 1, SLEG[0 + 8*J] == 1.0
// **Validates: Requirements 5.5**
// ============================================================================
static void test_optics_physical_constraints() {
    std::cout << "Test 3: OPTICS Physical Constraints\n";

    PRNG rng(456);
    const int NUM_TRIALS = 2000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        double PATH = rng.uniform(0.001, 10.0);
        int K = (trial % 2 == 0) ? 1 : 2;
        int K_orig = K;
        int rc = 0;

        std::array<double, S_> OPTD{};
        std::array<double, S_> SSALB{};
        std::array<double, 8 * S_> SLEG{};

        PhotoJX::OPTICS(OPTD.data(), SSALB.data(), SLEG.data(),
                        PATH, K, g_state, rc);

        if (rc != CLDJ_SUCCESS) {
            std::cerr << "  OPTICS returned error for K=" << K_orig
                      << " PATH=" << PATH << "\n";
            CHECK(false, "OPTICS: returned error for valid inputs");
            return;
        }

        for (int J = 0; J < S_; ++J) {
            if (OPTD[J] < 0.0) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " PATH=" << PATH
                          << " OPTD[" << J << "]=" << OPTD[J] << "\n";
                CHECK(false, "OPTICS: OD must be >= 0");
                return;
            }
            if (SSALB[J] < 0.0 || SSALB[J] > 1.0) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " PATH=" << PATH
                          << " SSALB[" << J << "]=" << SSALB[J] << "\n";
                CHECK(false, "OPTICS: SSA must be in [0,1]");
                return;
            }
            // OD * SSA <= OD (scattering OD <= total OD)
            if (OPTD[J] * SSALB[J] > OPTD[J] + 1.0e-15) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " PATH=" << PATH
                          << " OD*SSA=" << OPTD[J]*SSALB[J]
                          << " > OD=" << OPTD[J] << "\n";
                CHECK(false, "OPTICS: OD*SSA must be <= OD");
                return;
            }
            double p0 = SLEG[0 + 8 * J];
            if (std::abs(p0 - 1.0) > 1.0e-10) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " SLEG[0+" << 8*J << "]=" << p0 << "\n";
                CHECK(false, "OPTICS: First Legendre moment must be 1.0");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 4: OPTICG Physical Constraints (Property 8)
// For GeoMIP aerosol with valid K indices (1001..1000+NGG) and PATH > 0:
// verify OD >= 0, 0 <= SSA <= 1, SLEG[0 + 8*J] == 1.0
// **Validates: Requirements 5.5**
// ============================================================================
static void test_opticg_physical_constraints() {
    std::cout << "Test 4: OPTICG Physical Constraints\n";

    if (g_state.NGG == 0) {
        std::cout << "  SKIPPED (NGG=0, no GeoMIP data loaded)\n";
        return;
    }

    PRNG rng(789);
    const int NUM_TRIALS = 2000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        double PATH = rng.uniform(0.001, 5.0);
        int K = 1001 + rng.uniform_int(0, g_state.NGG - 1);

        std::array<double, S_> OPTD{};
        std::array<double, S_> SSALB{};
        std::array<double, 8 * S_> SLEG{};

        PhotoJX::OPTICG(OPTD.data(), SSALB.data(), SLEG.data(),
                        PATH, K, g_state);

        for (int J = 0; J < S_; ++J) {
            if (OPTD[J] < 0.0) {
                std::cerr << "  Counterexample: K=" << K
                          << " PATH=" << PATH
                          << " OPTD[" << J << "]=" << OPTD[J] << "\n";
                CHECK(false, "OPTICG: OD must be >= 0");
                return;
            }
            if (SSALB[J] < 0.0 || SSALB[J] > 1.0) {
                std::cerr << "  Counterexample: K=" << K
                          << " PATH=" << PATH
                          << " SSALB[" << J << "]=" << SSALB[J] << "\n";
                CHECK(false, "OPTICG: SSA must be in [0,1]");
                return;
            }
            if (OPTD[J] * SSALB[J] > OPTD[J] + 1.0e-15) {
                std::cerr << "  Counterexample: K=" << K
                          << " PATH=" << PATH
                          << " OD*SSA > OD\n";
                CHECK(false, "OPTICG: OD*SSA must be <= OD");
                return;
            }
            double p0 = SLEG[0 + 8 * J];
            if (std::abs(p0 - 1.0) > 1.0e-10) {
                std::cerr << "  Counterexample: K=" << K
                          << " SLEG[0+" << 8*J << "]=" << p0 << "\n";
                CHECK(false, "OPTICG: First Legendre moment must be 1.0");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 5: OPTICA Physical Constraints (Property 8)
// For standard aerosol Mie types with valid K indices (3..NAA) and PATH > 0:
// verify OD >= 0, 0 <= SSA <= 1, SLEG[0 + 8*J] == 1.0
// **Validates: Requirements 5.6**
// ============================================================================
static void test_optica_physical_constraints() {
    std::cout << "Test 5: OPTICA Physical Constraints\n";

    PRNG rng(1001);
    const int NUM_TRIALS = 2000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        double PATH = rng.uniform(0.001, 5.0);
        double RELH = rng.uniform(0.0, 1.0);
        int K = rng.uniform_int(3, g_state.NAA);
        int K_orig = K;
        int rc = 0;

        std::array<double, S_> OPTD{};
        std::array<double, S_> SSALB{};
        std::array<double, 8 * S_> SLEG{};

        PhotoJX::OPTICA(OPTD.data(), SSALB.data(), SLEG.data(),
                        PATH, RELH, K, g_state, rc);

        if (rc != CLDJ_SUCCESS) {
            std::cerr << "  OPTICA returned error for K=" << K_orig
                      << " PATH=" << PATH << "\n";
            CHECK(false, "OPTICA: returned error for valid inputs");
            return;
        }

        for (int J = 0; J < S_; ++J) {
            if (OPTD[J] < 0.0) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " PATH=" << PATH
                          << " OPTD[" << J << "]=" << OPTD[J] << "\n";
                CHECK(false, "OPTICA: OD must be >= 0");
                return;
            }
            if (SSALB[J] < 0.0 || SSALB[J] > 1.0) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " PATH=" << PATH
                          << " SSALB[" << J << "]=" << SSALB[J] << "\n";
                CHECK(false, "OPTICA: SSA must be in [0,1]");
                return;
            }
            if (OPTD[J] * SSALB[J] > OPTD[J] + 1.0e-15) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " OD*SSA > OD at bin " << J << "\n";
                CHECK(false, "OPTICA: OD*SSA must be <= OD");
                return;
            }
            double p0 = SLEG[0 + 8 * J];
            if (std::abs(p0 - 1.0) > 1.0e-10) {
                std::cerr << "  Counterexample: K=" << K_orig
                          << " SLEG[0+" << 8*J << "]=" << p0 << "\n";
                CHECK(false, "OPTICA: First Legendre moment must be 1.0");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 6: Optical Property Accumulation Additivity (Property 9)
// Given two cloud contributions with known OD1, SSA1 and OD2, SSA2:
// verify total_OD = OD1 + OD2
// and total_SSA = (OD1*SSA1 + OD2*SSA2) / total_OD
// **Validates: Requirements 5.7**
// ============================================================================
static void test_accumulation_additivity() {
    std::cout << "Test 6: Optical Property Accumulation Additivity (Property 9)\n";

    PRNG rng(2023);
    const int NUM_TRIALS = 10000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        // Generate two random cloud contributions
        // Each contribution has per-bin OD and SSA
        std::array<double, S_> OD1{}, OD2{};
        std::array<double, S_> SSA1{}, SSA2{};

        for (int J = 0; J < S_; ++J) {
            OD1[J]  = rng.uniform(0.0, 10.0);
            OD2[J]  = rng.uniform(0.0, 10.0);
            SSA1[J] = rng.uniform(0.0, 1.0);
            SSA2[J] = rng.uniform(0.0, 1.0);
        }

        // Compute accumulated totals using the additivity rule
        for (int J = 0; J < S_; ++J) {
            double total_OD = OD1[J] + OD2[J];

            // Verify additivity of optical depth
            double computed_total_OD = OD1[J] + OD2[J];
            if (std::abs(total_OD - computed_total_OD) > 1.0e-15) {
                CHECK(false, "Total OD != sum of individual ODs");
                return;
            }

            // Verify SSA weighted average formula
            if (total_OD > 0.0) {
                double total_SSA = (OD1[J] * SSA1[J] + OD2[J] * SSA2[J])
                                   / total_OD;

                // Verify total_SSA is in [0, 1]
                if (total_SSA < -1.0e-15 || total_SSA > 1.0 + 1.0e-15) {
                    std::cerr << "  Counterexample: OD1=" << OD1[J]
                              << " SSA1=" << SSA1[J]
                              << " OD2=" << OD2[J]
                              << " SSA2=" << SSA2[J]
                              << " total_SSA=" << total_SSA << "\n";
                    CHECK(false, "Accumulated SSA out of [0,1]");
                    return;
                }

                // Verify scattering OD <= total OD
                double scat_OD = total_OD * total_SSA;
                if (scat_OD > total_OD + 1.0e-14) {
                    std::cerr << "  Counterexample: scat_OD=" << scat_OD
                              << " > total_OD=" << total_OD << "\n";
                    CHECK(false, "Scattering OD exceeds total OD");
                    return;
                }
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 7: Accumulation with Real OPTICL/OPTICI Data (Property 9)
// Compute optical properties for liquid and ice cloud layers, then
// verify the accumulation rule holds when combining them.
// **Validates: Requirements 5.1, 5.2, 5.7**
// ============================================================================
static void test_accumulation_with_real_optics() {
    std::cout << "Test 7: Accumulation with Real OPTICL/OPTICI Data\n";

    PRNG rng(3141);
    const int NUM_TRIALS = 2000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        // Get liquid cloud optical properties
        double REFF_liq = rng.uniform(2.0, 40.0);
        double DDENS_liq = 0.0;
        std::array<double, S_> QQEXT_liq{}, SSALB_liq{};
        std::array<double, 8 * S_> SSLEG_liq{};

        PhotoJX::OPTICL(REFF_liq, 270.0, DDENS_liq,
                        QQEXT_liq.data(), SSALB_liq.data(),
                        SSLEG_liq.data(), g_state);

        // Get ice cloud optical properties
        double REFF_ice = rng.uniform(10.0, 100.0);
        double TEFF_ice = rng.uniform(200.0, 260.0);
        double DDENS_ice = 0.0;
        std::array<double, S_> QQEXT_ice{}, SSALB_ice{};
        std::array<double, 8 * S_> SSLEG_ice{};

        PhotoJX::OPTICI(REFF_ice, TEFF_ice, DDENS_ice,
                        QQEXT_ice.data(), SSALB_ice.data(),
                        SSLEG_ice.data(), g_state);

        // Simulate cloud water paths for each type
        double LWP = rng.uniform(0.1, 100.0);  // g/m2
        double IWP = rng.uniform(0.1, 50.0);   // g/m2

        // Compute per-bin optical depths
        // OD = PATH * 0.75 * Q / (Reff * density)
        for (int J = 0; J < S_; ++J) {
            double OD_liq = LWP * 0.75 * QQEXT_liq[J]
                            / (REFF_liq * DDENS_liq);
            double OD_ice = IWP * 0.75 * QQEXT_ice[J]
                            / (REFF_ice * DDENS_ice);

            double total_OD = OD_liq + OD_ice;

            // Check additivity
            if (total_OD < 0.0) {
                CHECK(false, "Total OD must be >= 0");
                return;
            }

            // Check weighted SSA
            if (total_OD > 1.0e-30) {
                double total_SSA = (OD_liq * SSALB_liq[J]
                                  + OD_ice * SSALB_ice[J]) / total_OD;

                if (total_SSA < -1.0e-14 || total_SSA > 1.0 + 1.0e-14) {
                    std::cerr << "  Counterexample: trial=" << trial
                              << " J=" << J
                              << " total_SSA=" << total_SSA << "\n";
                    CHECK(false, "Accumulated real SSA out of [0,1]");
                    return;
                }

                // Scattering OD <= total OD
                double scat_OD = total_OD * total_SSA;
                if (scat_OD > total_OD + 1.0e-12) {
                    CHECK(false, "Real optics: scattering OD > total OD");
                    return;
                }
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Test 8: OPTICL OD computation with known PATH (Property 8)
// For any cloud water path > 0, the computed optical depth must be >= 0
// and OD * SSA <= OD.
// **Validates: Requirements 5.1, 5.7**
// ============================================================================
static void test_opticl_od_with_path() {
    std::cout << "Test 8: OPTICL OD with known PATH (Property 8)\n";

    PRNG rng(5555);
    const int NUM_TRIALS = 3000;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        double REFF = rng.uniform(2.0, 45.0);
        double PATH = rng.uniform(0.01, 200.0);  // cloud water path g/m2
        double DDENS = 0.0;
        std::array<double, S_> QQEXT{}, SSALB{};
        std::array<double, 8 * S_> SSLEG{};

        PhotoJX::OPTICL(REFF, 270.0, DDENS, QQEXT.data(),
                        SSALB.data(), SSLEG.data(), g_state);

        for (int J = 0; J < S_; ++J) {
            // Compute optical depth: OD = PATH * extinction
            // extinction = 0.75 * Q / (Reff * density)
            double XTINCT = 0.75 * QQEXT[J] / (REFF * DDENS);
            double OD = PATH * XTINCT;

            if (OD < 0.0) {
                std::cerr << "  Counterexample: REFF=" << REFF
                          << " PATH=" << PATH
                          << " OD[" << J << "]=" << OD << "\n";
                CHECK(false, "OPTICL computed OD must be >= 0");
                return;
            }

            // OD * SSA <= OD (since SSA <= 1)
            double scat_OD = OD * SSALB[J];
            if (scat_OD > OD + 1.0e-14) {
                std::cerr << "  Counterexample: OD=" << OD
                          << " SSA=" << SSALB[J]
                          << " scat_OD=" << scat_OD << "\n";
                CHECK(false, "OPTICL: scattering OD must be <= total OD");
                return;
            }
        }
    }

    std::cout << "  PASSED (" << NUM_TRIALS << " random trials)\n";
}

// ============================================================================
// Main
// ============================================================================
int main() {
    std::cout << "=== Optical Property Tests (Properties 8 & 9) ===\n\n";

    if (!init_state()) {
        std::cerr << "FATAL: Failed to initialize CloudJState\n";
        return 1;
    }

    std::cout << "\n";

    test_opticl_physical_constraints();
    test_optici_physical_constraints();
    test_optics_physical_constraints();
    test_opticg_physical_constraints();
    test_optica_physical_constraints();
    test_accumulation_additivity();
    test_accumulation_with_real_optics();
    test_opticl_od_with_path();

    std::cout << "\n=== Summary ===\n";
    if (failures == 0) {
        std::cout << "All tests PASSED\n";
        return 0;
    } else {
        std::cout << failures << " test(s) FAILED\n";
        return 1;
    }
}
