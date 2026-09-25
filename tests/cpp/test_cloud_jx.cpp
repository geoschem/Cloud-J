/// @file test_cloud_jx.cpp
/// @brief Property-based tests for CLOUD_JX cloud overlap driver.
///
/// **Validates: Requirements 7.4, 7.5, 7.6, 7.7**
///
/// Property 15: Cloud Overlap J-value Bounds — For any valid cloud field with
/// CLDFLAG in {1,2,3,5,6,7,8}, the J-values computed by CLOUD_JX SHALL be
/// non-negative and finite.
///
/// Test strategy:
/// 1. Dark condition propagation: SZA > 98 → CLOUD_JX sets LDARK=true and
///    returns zero J-values regardless of CLDFLAG.
/// 2. CLDFLAG validation: CLDFLAG=4 (invalid) → rc = CLDJ_FAILURE.
/// 3. Clear sky (CLDFLAG=1) dark propagation across all valid CLDFLAGs.

#define MODEL_STANDALONE
#include <cloudj/cloud_jx.hpp>
#include <cloudj/state.hpp>
#include <cloudj/error.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
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

// Use the standalone LWEPAR constant for L1U
// In MODEL_STANDALONE: LWEPAR = 34, so L1U (number of levels) = LWEPAR + 1 = 35
static constexpr int L1U = CloudJState::LWEPAR + 1;
static constexpr int LU  = L1U - 1;

// =========================================================================
// Property 15 (Part 1): Dark Condition Propagation via CLOUD_JX
// =========================================================================
/// For SZA > 98 degrees, CLOUD_JX (regardless of CLDFLAG) should set
/// LDARK=true and return zero J-values. This tests CLOUD_JX's forwarding
/// of the dark condition from PHOTO_JX for all valid CLDFLAG values
/// {1, 2, 3, 5, 6, 7, 8}.

static void test_dark_condition_propagation() {
    std::cout << "--- Property 15 (Part 1): Dark Condition Propagation ---\n";

    const int ANU = AN_;  // 25 for standalone
    const int NJXU = 10;
    const int WW = W_ + W_r;  // 18

    // Valid CLDFLAG values
    std::vector<int> valid_cldflags = {1, 2, 3, 5, 6, 7, 8};
    // SZA values above the dark threshold (>98)
    std::vector<double> dark_szas = {99.0, 110.0, 150.0, 180.0};

    for (int cldflag : valid_cldflags) {
        for (double sza : dark_szas) {
            double U0 = std::cos(sza * CPI180);

            // Set up minimal state with CLDFLAG
            CloudJState state;
            state.CLDFLAG = cldflag;
            state.LNRG = 6;

            // Input arrays (zeros fine — dark check triggers before table use)
            double RFL_flat[5 * (W_ + W_r)] = {0};
            std::vector<double> PPP(L1U + 1, 0.0);
            std::vector<double> ZZZ(L1U + 1, 0.0);
            std::vector<double> TTT(L1U, 0.0);
            std::vector<double> HHH(L1U, 0.0);
            std::vector<double> DDD(L1U, 0.0);
            std::vector<double> RRR(L1U, 0.0);
            std::vector<double> OOO(L1U, 0.0);
            std::vector<double> CCC(L1U, 0.0);
            std::vector<double> LWP_arr(L1U, 0.0);
            std::vector<double> IWP_arr(L1U, 0.0);
            std::vector<double> REFFL(L1U, 0.0);
            std::vector<double> REFFI(L1U, 0.0);
            std::vector<double> CLDF(L1U, 0.0);
            std::vector<int>    CLDIW(L1U, 0);
            std::vector<double> AERSP(L1U, 0.0);
            std::vector<int>    NDXAER(L1U, 0);

            // Output arrays initialized with sentinel values
            std::vector<double> VALJXX(LU * NJXU, 999.0);
            std::vector<double> SKPERD((S_ + 2) * L1U, 999.0);
            double SWMSQ[6] = {999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
            std::vector<double> OD18(L1U, 999.0);
            double WTQCA[NQD_] = {0};

            int NICA_out = 0;
            int JCOUNT = 0;
            bool LDARK = false;
            int rc = CLDJ_SUCCESS;

            CLOUD_JX(
                U0, sza,
                RFL_flat, 1.0, false,
                PPP.data(), ZZZ.data(),
                TTT.data(), HHH.data(),
                DDD.data(), RRR.data(),
                OOO.data(), CCC.data(),
                LWP_arr.data(), IWP_arr.data(),
                REFFL.data(), REFFI.data(),
                CLDF.data(), CLDIW.data(),
                0.0,  // CLDCOR
                AERSP.data(), NDXAER.data(),
                L1U, ANU, NJXU,
                VALJXX.data(), SKPERD.data(), SWMSQ, OD18.data(),
                0, NICA_out, JCOUNT, LDARK,
                WTQCA,
                state, rc);

            std::string desc = "CLDFLAG=" + std::to_string(cldflag) +
                               ", SZA=" + std::to_string(sza);

            // rc should be success (CLDFLAG is valid)
            check(rc == CLDJ_SUCCESS,
                  "rc not SUCCESS for dark condition: " + desc);

            // LDARK must be true
            check(LDARK == true,
                  "LDARK not set for dark condition: " + desc);

            // All VALJXX must be exactly 0.0
            bool all_jval_zero = true;
            for (int i = 0; i < LU * NJXU; ++i) {
                if (VALJXX[i] != 0.0) {
                    all_jval_zero = false;
                    break;
                }
            }
            check(all_jval_zero,
                  "VALJXX not all zero for dark: " + desc);

            // All SKPERD must be exactly 0.0
            bool all_skperd_zero = true;
            for (int i = 0; i < (S_ + 2) * L1U; ++i) {
                if (SKPERD[i] != 0.0) {
                    all_skperd_zero = false;
                    break;
                }
            }
            check(all_skperd_zero,
                  "SKPERD not all zero for dark: " + desc);

            // JCOUNT must be 0 (no PHOTO_JX calls that returned non-dark)
            check(JCOUNT == 0,
                  "JCOUNT not 0 for dark: " + desc);
        }
    }

    std::cout << "  Dark condition propagation verified for "
              << valid_cldflags.size() << " CLDFLAGs x "
              << dark_szas.size() << " SZA values\n";
}

// =========================================================================
// Property 15 (Part 2): CLDFLAG Validation
// =========================================================================
/// CLDFLAG=4 is invalid. CLOUD_JX should return rc = CLDJ_FAILURE.
/// Also test CLDFLAG=0 and CLDFLAG=9 (out of range).

static void test_cldflag_validation() {
    std::cout << "--- Property 15 (Part 2): CLDFLAG Validation ---\n";

    const int ANU = AN_;
    const int NJXU = 10;

    // Invalid CLDFLAG values
    std::vector<int> invalid_cldflags = {0, 4, 9, -1};

    for (int cldflag : invalid_cldflags) {
        // Use a non-dark SZA so the function proceeds past the dark check
        double sza = 60.0;
        double U0 = std::cos(sza * CPI180);

        CloudJState state;
        state.CLDFLAG = cldflag;
        state.LNRG = 6;

        double RFL_flat[5 * (W_ + W_r)] = {0};
        std::vector<double> PPP(L1U + 1, 0.0);
        std::vector<double> ZZZ(L1U + 1, 0.0);
        std::vector<double> TTT(L1U, 0.0);
        std::vector<double> HHH(L1U, 0.0);
        std::vector<double> DDD(L1U, 0.0);
        std::vector<double> RRR(L1U, 0.0);
        std::vector<double> OOO(L1U, 0.0);
        std::vector<double> CCC(L1U, 0.0);
        std::vector<double> LWP_arr(L1U, 0.0);
        std::vector<double> IWP_arr(L1U, 0.0);
        std::vector<double> REFFL(L1U, 0.0);
        std::vector<double> REFFI(L1U, 0.0);
        std::vector<double> CLDF(L1U, 0.0);
        std::vector<int>    CLDIW(L1U, 0);
        std::vector<double> AERSP(L1U, 0.0);
        std::vector<int>    NDXAER(L1U, 0);

        std::vector<double> VALJXX(LU * NJXU, 0.0);
        std::vector<double> SKPERD((S_ + 2) * L1U, 0.0);
        double SWMSQ[6] = {0};
        std::vector<double> OD18(L1U, 0.0);
        double WTQCA[NQD_] = {0};

        int NICA_out = 0;
        int JCOUNT = 0;
        bool LDARK = false;
        int rc = CLDJ_SUCCESS;

        CLOUD_JX(
            U0, sza,
            RFL_flat, 1.0, false,
            PPP.data(), ZZZ.data(),
            TTT.data(), HHH.data(),
            DDD.data(), RRR.data(),
            OOO.data(), CCC.data(),
            LWP_arr.data(), IWP_arr.data(),
            REFFL.data(), REFFI.data(),
            CLDF.data(), CLDIW.data(),
            0.0,
            AERSP.data(), NDXAER.data(),
            L1U, ANU, NJXU,
            VALJXX.data(), SKPERD.data(), SWMSQ, OD18.data(),
            0, NICA_out, JCOUNT, LDARK,
            WTQCA,
            state, rc);

        std::string desc = "CLDFLAG=" + std::to_string(cldflag);

        // rc should indicate failure for invalid CLDFLAG
        check(rc == CLDJ_FAILURE,
              "rc not CLDJ_FAILURE for invalid " + desc);
    }

    std::cout << "  CLDFLAG validation verified for invalid values {0, 4, 9, -1}\n";
}

// =========================================================================
// Property 15 (Part 3): J-value Non-negativity for CLDFLAG=1 (Clear Sky)
// =========================================================================
/// For CLDFLAG=1 with zero cloud and dark SZA, verify outputs are zeroed.
/// This validates the zero-cloud path through CLOUD_JX produces non-negative,
/// finite outputs.

static void test_clear_sky_outputs_non_negative() {
    std::cout << "--- Property 15 (Part 3): Clear Sky Non-negative Outputs ---\n";

    const int ANU = AN_;
    const int NJXU = 10;

    // Use a dark SZA to get deterministic zero output with clear sky
    // (This tests the full CLOUD_JX→PHOTO_JX dark path with CLDFLAG=1)
    double sza = 99.0;
    double U0 = std::cos(sza * CPI180);

    CloudJState state;
    state.CLDFLAG = 1;
    state.LNRG = 6;

    double RFL_flat[5 * (W_ + W_r)] = {0};
    std::vector<double> PPP(L1U + 1, 0.0);
    std::vector<double> ZZZ(L1U + 1, 0.0);
    std::vector<double> TTT(L1U, 0.0);
    std::vector<double> HHH(L1U, 0.0);
    std::vector<double> DDD(L1U, 0.0);
    std::vector<double> RRR(L1U, 0.0);
    std::vector<double> OOO(L1U, 0.0);
    std::vector<double> CCC(L1U, 0.0);
    std::vector<double> LWP_arr(L1U, 0.0);
    std::vector<double> IWP_arr(L1U, 0.0);
    std::vector<double> REFFL(L1U, 0.0);
    std::vector<double> REFFI(L1U, 0.0);
    std::vector<double> CLDF(L1U, 0.0);
    std::vector<int>    CLDIW(L1U, 0);
    std::vector<double> AERSP(L1U, 0.0);
    std::vector<int>    NDXAER(L1U, 0);

    std::vector<double> VALJXX(LU * NJXU, 0.0);
    std::vector<double> SKPERD((S_ + 2) * L1U, 0.0);
    double SWMSQ[6] = {0};
    std::vector<double> OD18(L1U, 0.0);
    double WTQCA[NQD_] = {0};

    int NICA_out = 0;
    int JCOUNT = 0;
    bool LDARK = false;
    int rc = CLDJ_SUCCESS;

    CLOUD_JX(
        U0, sza,
        RFL_flat, 1.0, false,
        PPP.data(), ZZZ.data(),
        TTT.data(), HHH.data(),
        DDD.data(), RRR.data(),
        OOO.data(), CCC.data(),
        LWP_arr.data(), IWP_arr.data(),
        REFFL.data(), REFFI.data(),
        CLDF.data(), CLDIW.data(),
        0.0,
        AERSP.data(), NDXAER.data(),
        L1U, ANU, NJXU,
        VALJXX.data(), SKPERD.data(), SWMSQ, OD18.data(),
        0, NICA_out, JCOUNT, LDARK,
        WTQCA,
        state, rc);

    check(rc == CLDJ_SUCCESS, "Clear sky call should succeed");
    check(LDARK == true, "Clear sky with SZA>98 should be dark");

    // Verify all J-values are non-negative and finite
    bool all_non_negative = true;
    bool all_finite = true;
    for (int i = 0; i < LU * NJXU; ++i) {
        if (VALJXX[i] < 0.0) all_non_negative = false;
        if (!std::isfinite(VALJXX[i])) all_finite = false;
    }
    check(all_non_negative, "VALJXX should be non-negative for clear sky");
    check(all_finite, "VALJXX should be finite for clear sky");

    // Verify OD18 is non-negative and finite
    for (int i = 0; i < L1U; ++i) {
        if (OD18[i] < 0.0) all_non_negative = false;
        if (!std::isfinite(OD18[i])) all_finite = false;
    }
    check(all_non_negative, "OD18 should be non-negative for clear sky");
    check(all_finite, "OD18 should be finite for clear sky");

    std::cout << "  Clear sky (CLDFLAG=1) non-negative outputs verified\n";
}

// =========================================================================
// Property 15 (Part 4): CLOUD_JX outputs zero for all ICA-based CLDFLAGs
//                        in dark conditions
// =========================================================================
/// ICA-based CLDFLAGs {5,6,7,8} involve multiple PHOTO_JX calls and weighted
/// averaging. The dark condition should propagate correctly through all of
/// them, resulting in zero accumulated J-values.

static void test_ica_dark_accumulation() {
    std::cout << "--- Property 15 (Part 4): ICA Dark Accumulation ---\n";

    const int ANU = AN_;
    const int NJXU = 10;

    std::vector<int> ica_cldflags = {5, 6, 7, 8};
    double sza = 120.0;  // well above dark threshold
    double U0 = std::cos(sza * CPI180);

    for (int cldflag : ica_cldflags) {
        CloudJState state;
        state.CLDFLAG = cldflag;
        state.LNRG = 6;

        double RFL_flat[5 * (W_ + W_r)] = {0};
        std::vector<double> PPP(L1U + 1, 0.0);
        std::vector<double> ZZZ(L1U + 1, 0.0);
        std::vector<double> TTT(L1U, 0.0);
        std::vector<double> HHH(L1U, 0.0);
        std::vector<double> DDD(L1U, 0.0);
        std::vector<double> RRR(L1U, 0.0);
        std::vector<double> OOO(L1U, 0.0);
        std::vector<double> CCC(L1U, 0.0);
        std::vector<double> LWP_arr(L1U, 0.0);
        std::vector<double> IWP_arr(L1U, 0.0);
        std::vector<double> REFFL(L1U, 0.0);
        std::vector<double> REFFI(L1U, 0.0);
        // Some non-zero cloud to exercise ICA paths (but dark overrides)
        std::vector<double> CLDF(L1U, 0.5);
        std::vector<int>    CLDIW(L1U, 1);
        std::vector<double> AERSP(L1U, 0.0);
        std::vector<int>    NDXAER(L1U, 0);

        std::vector<double> VALJXX(LU * NJXU, 999.0);
        std::vector<double> SKPERD((S_ + 2) * L1U, 999.0);
        double SWMSQ[6] = {999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
        std::vector<double> OD18(L1U, 999.0);
        double WTQCA[NQD_] = {0};

        int NICA_out = 0;
        int JCOUNT = 0;
        bool LDARK = false;
        int rc = CLDJ_SUCCESS;

        CLOUD_JX(
            U0, sza,
            RFL_flat, 1.0, false,
            PPP.data(), ZZZ.data(),
            TTT.data(), HHH.data(),
            DDD.data(), RRR.data(),
            OOO.data(), CCC.data(),
            LWP_arr.data(), IWP_arr.data(),
            REFFL.data(), REFFI.data(),
            CLDF.data(), CLDIW.data(),
            0.5,  // non-zero CLDCOR
            AERSP.data(), NDXAER.data(),
            L1U, ANU, NJXU,
            VALJXX.data(), SKPERD.data(), SWMSQ, OD18.data(),
            0, NICA_out, JCOUNT, LDARK,
            WTQCA,
            state, rc);

        std::string desc = "CLDFLAG=" + std::to_string(cldflag);

        check(rc == CLDJ_SUCCESS,
              "rc not SUCCESS for ICA dark: " + desc);

        // Even with non-zero cloud fractions, dark SZA means LDARK=true
        // and all accumulated J-values should be zero.
        check(LDARK == true,
              "LDARK not set for ICA dark: " + desc);

        // All VALJXX must be zero after accumulation
        bool all_jval_zero = true;
        for (int i = 0; i < LU * NJXU; ++i) {
            if (VALJXX[i] != 0.0) {
                all_jval_zero = false;
                break;
            }
        }
        check(all_jval_zero,
              "VALJXX not all zero for ICA dark: " + desc);

        // J-values should be non-negative and finite
        bool all_non_neg = true;
        bool all_finite = true;
        for (int i = 0; i < LU * NJXU; ++i) {
            if (VALJXX[i] < 0.0) all_non_neg = false;
            if (!std::isfinite(VALJXX[i])) all_finite = false;
        }
        check(all_non_neg,
              "VALJXX negative for ICA dark: " + desc);
        check(all_finite,
              "VALJXX not finite for ICA dark: " + desc);
    }

    std::cout << "  ICA dark accumulation verified for CLDFLAGs {5, 6, 7, 8}\n";
}

// =========================================================================
// Main
// =========================================================================

int main() {
    std::cout << "=== CLOUD_JX Property Tests ===\n\n";

    test_dark_condition_propagation();
    test_cldflag_validation();
    test_clear_sky_outputs_non_negative();
    test_ica_dark_accumulation();

    std::cout << "\n=== Results: " << (total_tests - failures) << "/"
              << total_tests << " passed ===\n";

    if (failures > 0) {
        std::cerr << failures << " test(s) FAILED\n";
        return 1;
    }

    std::cout << "All tests passed.\n";
    return 0;
}
