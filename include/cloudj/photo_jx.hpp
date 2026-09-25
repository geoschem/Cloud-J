#ifndef CLOUDJ_PHOTO_JX_HPP
#define CLOUDJ_PHOTO_JX_HPP

#include <cloudj/state.hpp>
#include <cloudj/error.hpp>
#include <cloudj/cross_sections.hpp>
#include <cloudj/fast_math.hpp>
#include <cloudj/radiative_solver.hpp>
#include <cloudj/photolysis.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <experimental/mdspan.hpp>

namespace CloudJ {
namespace PhotoJX {

// Convert a double to Fortran's E9.2 scientific editing (width 9, e.g.
// " 4.71E-05"). Used by the LPRTJ diagnostic print in PHOTO_JX.
inline std::string format_fortran_e9_2(double val) {
    if (std::abs(val) < 1e-99) return " 0.00E+00";
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%9.2E", val);
    std::string s(buf);
    // std::printf may emit a 3-digit exponent (e.g. 1.00E-005 on some
    // platforms); Fortran E9.2 always uses a sign plus exactly 2 digits.
    const size_t epos = s.find('E');
    if (epos != std::string::npos) {
        const std::string mant = s.substr(0, epos);
        const std::string exp = s.substr(epos + 1);
        char fixed[32];
        if (exp[0] == '-') {
            std::snprintf(fixed, sizeof(fixed), "%sE-%02d", mant.c_str(),
                          std::atoi(exp.c_str() + 1));
        } else {
            std::snprintf(fixed, sizeof(fixed), "%sE+%02d", mant.c_str(),
                          std::atoi(exp.c_str() + 1));
        }
        s = fixed;
        while (s.length() < 9) s = " " + s;
    }
    return s;
}

// =========================================================================
// SPHERE1F - Flat Earth Air Mass Factors
// =========================================================================
/// Compute flat-Earth air mass factors: AMF = 1/cos(SZA) for all layers.
/// Needed for testing flat-disk errors against spherical geometry.
///
/// Ported from SPHERE1F in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param U0   cos(solar zenith angle)
/// @param RAD  radius of Earth mean sea level (cm) [unused in flat-Earth]
/// @param ZHL  height (cm) of bottom edge of CTM levels, dimension L1U+1
/// @param ZZHT scale height (cm) used above top of CTM [unused in flat-Earth]
/// @param AMF  output air mass factor array, dimension (L1U+1)*(L1U+1),
///             stored column-major: AMF[J + (L1U+1)*L] = AMF(J,L)
/// @param L1U  number of CTM levels + 1 (dimension of CTM)
inline void SPHERE1F(double U0, double RAD,
                     const double* ZHL, double ZZHT,
                     double* AMF, int L1U)
{
    (void)RAD;   // unused in flat-Earth
    (void)ZHL;   // unused in flat-Earth
    (void)ZZHT;  // unused in flat-Earth

    const int LTOP = L1U;
    const int dim = L1U + 1;  // array dimension (Fortran 1:L1U+1)

    // Initialize AMF to zero
    for (int i = 0; i < dim * dim; ++i) {
        AMF[i] = 0.0;
    }

    // AMF(LTOP+1, LTOP+1) = 1.0 (Fortran indexing: top-of-atmos marker)
    // In C zero-based: AMF[LTOP + dim*LTOP] = 1.0
    AMF[LTOP + dim * LTOP] = 1.0;

    if (U0 > 0.0) {
        double PATH0 = 1.0 / U0;

        // Fortran: do L=1,LTOP; do J=L,LTOP; AMF(J,L)=PATH0
        // C zero-based: L -> L-1, J -> J-1
        // AMF(J,L) in Fortran = AMF[(J-1) + dim*(L-1)] in C
        for (int L = 0; L < LTOP; ++L) {
            for (int J = L; J < LTOP; ++J) {
                AMF[J + dim * L] = PATH0;
            }
        }
    }
}

// =========================================================================
// SPHERE1N - Spherical Straight-Line Path Air Mass Factors
// =========================================================================
/// Compute spherical geometry air mass factors using straight-line paths
/// (no atmospheric refraction).
///
/// Ported from SPHERE1N in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param U0   cos(solar zenith angle)
/// @param RAD  radius of Earth mean sea level (cm)
/// @param ZHL  height (cm) of bottom edge of CTM levels, dimension L1U+1
/// @param ZZHT scale height (cm) used above top of CTM
/// @param AMF  output air mass factor array, dimension (L1U+1)*(L1U+1),
///             stored column-major: AMF[J + (L1U+1)*L] = AMF(J,L)
/// @param L1U  number of CTM levels + 1 (dimension of CTM)
inline void SPHERE1N(double U0, double RAD,
                     const double* ZHL, double ZZHT,
                     double* AMF, int L1U)
{
    const int LTOP = L1U;
    const int dim = L1U + 1;  // array dimension

    // Local arrays (dynamically sized for generality)
    // Using VLAs via vectors for stack safety with large L1U
    double RZ[128], DIVZ[128], RATZ[128];  // max supported levels
    // If L1U+1 > 128, this would need dynamic allocation.
    // Cloud-J typically has L1U <= 58, so 128 is safe.

    // Initialize AMF to zero
    for (int i = 0; i < dim * dim; ++i) {
        AMF[i] = 0.0;
    }

    // Compute RZ (radius at each level edge)
    // Fortran: RZ(L) = RAD + ZHL(L) for L=1..LTOP
    // C zero-based: RZ[L] = RAD + ZHL[L] for L=0..LTOP-1
    for (int L = 0; L < LTOP; ++L) {
        RZ[L] = RAD + ZHL[L];
    }
    RZ[LTOP] = RZ[LTOP - 1] + ZZHT;

    // Compute DIVZ and RATZ
    for (int L = 0; L < LTOP; ++L) {
        DIVZ[L] = 1.0 / (RZ[L + 1] - RZ[L]);
        RATZ[L] = RZ[L] / RZ[L + 1];
    }

    double CA0 = U0;
    double A0 = std::acos(CA0);
    double SA0 = std::sin(A0);
    double R0 = RZ[0];

    // AMF(LTOP+1, LTOP+1) = 1.0
    AMF[LTOP + dim * LTOP] = 1.0;

    if (CA0 >= 0.0) {
        // Surface in direct sunlight
        // For each level L (Fortran L=1..LTOP, C L=0..LTOP-1):
        //   trace ray upward from level L through layers J=L..LTOP
        for (int L = 0; L < LTOP; ++L) {
            double SA1 = SA0;  // starting angle at radius point L
            double A1 = std::asin(SA1);
            double CA1 = std::cos(A1);

            for (int J = L; J < LTOP; ++J) {
                double SA2 = SA1 * RATZ[J];  // RZ[J]/RZ[J+1]
                double A2 = std::asin(SA2);
                double CA2 = std::cos(A2);
                double PATH = RZ[J + 1] * CA2 - RZ[J] * CA1;
                AMF[J + dim * L] = PATH * DIVZ[J];
                SA1 = SA2;
                CA1 = CA2;
                A1 = A2;
            }
        }
    } else {
        // Surface dark: search upward to find sunlit levels
        // Fortran: do L=2,LTOP+1 → C: L=1..LTOP
        for (int L = 1; L <= LTOP; ++L) {
            if (SA0 * RZ[L] > R0) {
                // This level is above the terminator
                double SA1 = SA0;
                double CA1 = CA0;
                int JUP = 0;
                bool crossed_terminator = false;

                // Trace ray downward from level L to find terminator
                // Fortran: do J=L-1,1,-1 → C: J=L-1..0 (layer indices)
                for (int J = L - 1; J >= 0; --J) {
                    if (SA1 * RZ[J + 1] < RZ[J]) {
                        // Path from R(J+1) down to R(J)
                        double SA2 = SA1 / RATZ[J];  // SA1 * RZ[J+1]/RZ[J]
                        double A2 = std::asin(SA2);
                        double CA2 = -std::cos(A2);
                        double PATH = RZ[J + 1] * CA2 - RZ[J] * CA1;
                        AMF[J + dim * L] += PATH * DIVZ[J];
                        SA1 = SA2;
                        CA1 = CA2;
                        (void)A2;  // angle carried by SA1/CA1
                    } else {
                        // Path across terminator (CA=0) in layer J
                        double PATH = -2.0 * CA1 * RZ[J + 1];
                        AMF[J + dim * L] += PATH * DIVZ[J];
                        CA1 = -CA1;
                        JUP = J + 1;
                        crossed_terminator = true;
                        break;
                    }
                }

                if (!crossed_terminator) {
                    // Ray went all the way down without crossing terminator
                    // (shouldn't normally happen for dark surface)
                    continue;
                }

                // Start back up from JUP
                // Fortran: do J=JUP,LTOP → C: J=JUP..LTOP-1
                for (int J = JUP; J < LTOP; ++J) {
                    double SA2 = SA1 * RATZ[J];
                    double A2 = std::asin(SA2);
                    double CA2 = std::cos(A2);
                    double PATH = RZ[J + 1] * CA2 - RZ[J] * CA1;
                    AMF[J + dim * L] += PATH * DIVZ[J];
                    SA1 = SA2;
                    CA1 = CA2;
                    (void)A2;  // angle carried by SA1/CA1
                }
            }
            // else: level still in Earth's shadow, AMF stays 0
        }
    }
}

// =========================================================================
// SPHERE1R - Spherical Air Mass Factors with Atmospheric Refraction
// =========================================================================
/// Compute spherical geometry air mass factors with atmospheric refraction
/// correction using Snell's law through atmospheric layers.
///
/// Ported from SPHERE1R in src/Core/cldj_fjx_sub_mod.F90.
///
/// Uses a density scale height of 8 km and refractive index of 1.000300
/// at sea level, scaling with atmospheric density.
///
/// @param U0   cos(solar zenith angle)
/// @param RAD  radius of Earth mean sea level (cm)
/// @param ZHL  height (cm) of bottom edge of CTM levels, dimension L1U+1
/// @param ZZHT scale height (cm) used above top of CTM
/// @param AMF  output air mass factor array, dimension (L1U+1)*(L1U+1),
///             stored column-major: AMF[J + (L1U+1)*L] = AMF(J,L)
/// @param L1U  number of CTM levels + 1 (dimension of CTM)
inline void SPHERE1R(double U0, double RAD,
                     const double* ZHL, double ZZHT,
                     double* AMF, int L1U)
{
    const int LTOP = L1U;
    const int dim = L1U + 1;

    // Local arrays
    double RZ[128], DIVZ[128], RATZ[128], RD[128], RN[128];
    double PATH1[128], PATH2[128], ZANG1[128];
    // 2D arrays for refracted tangent paths
    // ZANG(L,K) and ZAMF(L,K) - stored as ZANG[L*dim + K]
    double ZANG[128 * 128], ZAMF_arr[128 * 128];

    // Refraction parameters
    const double DDHT = 8.0e5;        // density scale height (cm) = 8 km
    const double REF0 = 300.0e-6;     // refractive index increment at sea level
    const double C90 = 1.570796326794897;  // pi/2

    // Initialize AMF to zero
    for (int i = 0; i < dim * dim; ++i) {
        AMF[i] = 0.0;
    }

    // Initialize 2D local arrays
    for (int i = 0; i < dim * dim; ++i) {
        ZANG[i] = 0.0;
        ZAMF_arr[i] = 0.0;
    }

    // Compute RZ, DIVZ, RATZ, RD, RN
    for (int L = 0; L < LTOP; ++L) {
        RZ[L] = RAD + ZHL[L];
    }
    RZ[LTOP] = RZ[LTOP - 1] + ZZHT;

    for (int L = 0; L < LTOP; ++L) {
        DIVZ[L] = 1.0 / (RZ[L + 1] - RZ[L]);
        RATZ[L] = RZ[L] / RZ[L + 1];
        RD[L] = RadiativeSolver::exp_eval(-(RZ[L] - RAD) / DDHT);
        RN[L] = 1.0 + REF0 * RD[L];
    }
    RD[LTOP] = 0.0;
    RN[LTOP] = 1.0;

    double CZA0 = U0;
    double ZA0 = std::acos(CZA0);
    double SZA0 = std::sin(ZA0);

    // AMF(LTOP+1, LTOP+1) = 1.0
    AMF[LTOP + dim * LTOP] = 1.0;

    if (U0 >= 0.0) {
        // ===== Surface in direct sunlight (ZA0 <= 90 deg) =====
        // Do first downward integration with refraction to get elevation angle,
        // then redo with corrected angle for each layer edge.
        for (int L = 0; L < LTOP; ++L) {
            double SRN0 = SZA0 * RZ[L] * RN[L];  // Snell invariant

            // First pass: compute elevation angle correction
            double SA0_top = SRN0 / (RZ[LTOP] * RN[LTOP]);
            ZANG1[LTOP] = std::asin(SA0_top);
            for (int i = 0; i < dim; ++i) PATH1[i] = 0.0;

            for (int K = LTOP - 1; K >= L; --K) {
                double SA1 = SRN0 / (RZ[K] * RN[K]);
                double SA2 = SA1 * RATZ[K];
                double A1 = std::asin(SA1);
                double A2 = std::asin(SA2);
                ZANG1[K] = ZANG1[K + 1] + A1 - A2;
            }

            // Correct zenith angle at lower edge L
            double ZA1 = ZA0 - (ZANG1[L] - ZA0);

            // Second pass with corrected angle: compute actual paths
            double SZA1 = std::sin(ZA1);
            double SRN1 = SZA1 * RZ[L] * RN[L];
            for (int i = 0; i < dim; ++i) PATH1[i] = 0.0;

            for (int K = LTOP - 1; K >= L; --K) {
                double SA1 = SRN1 / (RZ[K] * RN[K]);
                double SA2 = SA1 * RATZ[K];
                double A1 = std::asin(SA1);
                double A2 = std::asin(SA2);
                double CA1 = std::cos(A1);
                double CA2 = std::cos(A2);
                PATH1[K] = RZ[K + 1] * CA2 - RZ[K] * CA1;
            }

            // Store AMF
            for (int K = 0; K < LTOP; ++K) {
                AMF[K + dim * L] = PATH1[K] * DIVZ[K];
            }
        }
    } else {
        // ===== Surface dark (ZA0 > 90 deg) =====
        // Integrate refracted paths tangent at each radius RZ[L]

        // Build ZAMF and ZANG tables for tangent paths at each level
        for (int L = 0; L <= LTOP; ++L) {
            double SRN0 = RZ[L] * RN[L];  // Snell invariant for tangent at L
            double SA0_top = SRN0 / (RZ[LTOP] * RN[LTOP]);
            ZANG1[LTOP] = std::asin(SA0_top);
            for (int i = 0; i < dim; ++i) PATH1[i] = 0.0;

            // Trace downward to tangent layer L
            for (int K = LTOP - 1; K >= L; --K) {
                double SA1 = SRN0 / (RZ[K] * RN[K]);
                double SA2 = SA1 * RATZ[K];
                double A1 = std::asin(SA1);
                double A2 = std::asin(SA2);
                double CA1 = std::cos(A1);
                double CA2 = std::cos(A2);
                PATH1[K] = RZ[K + 1] * CA2 - RZ[K] * CA1;
                ZANG1[K] = ZANG1[K + 1] + A1 - A2;
            }

            // Store symmetric path and angle
            // ZAMF(L,K) = PATH1(K)*DIVZ(K), ZANG(L,K) = symmetric angle
            for (int K = 0; K < LTOP; ++K) {
                ZAMF_arr[L * dim + K] = PATH1[K] * DIVZ[K];
                ZANG[L * dim + K] = ZANG1[L] + ZANG1[L] - ZANG1[K];
            }
            int K = LTOP;
            ZANG[L * dim + K] = ZANG1[L] + ZANG1[L] - ZANG1[K];
        }

        // Now compute AMF for each level L0
        for (int L0 = 0; L0 <= LTOP; ++L0) {
            if (ZA0 < ZANG[L0 * dim + L0]) {
                // Pre-terminator: direct path exists
                if (L0 < LTOP) {
                    // Correct zenith angle using elevation angle
                    double ZA1 = ZA0 - (ZANG[L0 * dim + L0] - C90);
                    double SZA1 = std::sin(ZA1);
                    double SRN1 = SZA1 * RZ[L0] * RN[L0];
                    for (int i = 0; i < dim; ++i) PATH1[i] = 0.0;

                    for (int K = LTOP - 1; K >= L0; --K) {
                        double SA1 = SRN1 / (RZ[K] * RN[K]);
                        double SA2 = SA1 * RATZ[K];
                        double A1 = std::asin(SA1);
                        double A2 = std::asin(SA2);
                        double CA1_k = std::cos(A1);
                        double CA2_k = std::cos(A2);
                        PATH1[K] = RZ[K + 1] * CA2_k - RZ[K] * CA1_k;
                    }

                    for (int K = 0; K < LTOP; ++K) {
                        AMF[K + dim * L0] = PATH1[K] * DIVZ[K];
                    }
                }
            } else {
                // Post-terminator: ZA0 > ZANG(L0,L0)
                // Find interpolation bracket K0
                int K0 = -1;
                for (int K = 0; K < L0; ++K) {
                    if (ZA0 <= ZANG[K * dim + L0]) {
                        K0 = K;
                    }
                }

                if (K0 >= 0) {
                    // Interpolate between tangent paths K0 and K0+1
                    double F0 = (ZA0 - ZANG[(K0 + 1) * dim + L0]) /
                                (ZANG[K0 * dim + L0] - ZANG[(K0 + 1) * dim + L0]);

                    for (int K = 0; K < LTOP; ++K) {
                        PATH2[K] = F0 * ZAMF_arr[K0 * dim + K] +
                                   (1.0 - F0) * ZAMF_arr[(K0 + 1) * dim + K];
                    }

                    // Load pre-terminator paths
                    for (int K = 0; K < LTOP; ++K) {
                        AMF[K + dim * L0] = PATH2[K];
                    }

                    // Add post-terminator paths below level L0
                    for (int K = 0; K < L0; ++K) {
                        AMF[K + dim * L0] += PATH2[K];
                    }
                }
            }
        }
    }
}

// =========================================================================
// EXTRAL1 — Cloud Layer Insertion
// =========================================================================
/// Determines the number of extra sub-layers to insert at each level based on
/// cloud optical depth exceeding the ATAU0 threshold. Thick cloud layers are
/// subdivided using geometric factor ATAU between successive inserted layers.
///
/// Ported from EXTRAL1 in src/Core/cldj_fjx_sub_mod.F90.
///
/// The key parameters are:
///   ATAU  = factor increase from one layer to the next
///   ATAU0 = delta-TAU cut-off for cloud OD to insert a layer
///
/// For each layer (from top down), if its OD > ATAU0X:
///   JX = round(ln(1 + (ATAU-1)*DTAU/ATAU0X) / ln(ATAU))
///   JXTRA[L] = JX
///   ATAU0X *= ATAU^JX  (threshold grows geometrically)
///
/// Overflow check: if total expanded levels * 2 > NX, zero remaining JXTRA.
///
/// @param DTAU600  Optical depth per layer (generally 600nm, cloud+aerosol) [L1X]
/// @param L1X      Number of layers
/// @param NX       Max Mie scattering array size (N_=601)
/// @param ATAU     Geometric factor between successive inserted layers (default 1.05)
/// @param ATAU0    Minimum OD threshold for uppermost inserted layer (default 0.005)
/// @param JXTRA    Output: number of extra sub-layers at each level [L1X]
inline void EXTRAL1(
    const double* DTAU600,  // optical depth per layer [L1X]
    int L1X,                // number of layers
    int NX,                 // max Mie scattering array size (N_=601)
    double ATAU,            // geometric factor (default 1.05)
    double ATAU0,           // minimum OD threshold (default 0.005)
    int* JXTRA)             // output: extra layers per level [L1X]
{
    const double ATAULN = std::log(ATAU);
    double ATAU0X = ATAU0;

    // Loop from top (L1X-1) down to 0 — mirrors Fortran L=L1X,1,-1
    // Fortran 1-based indexing maps to C++ 0-based: DTAU600(L) → DTAU600[L-1]
    for (int L = L1X - 1; L >= 0; --L) {
        JXTRA[L] = 0;
        if (DTAU600[L] > ATAU0X) {
            double AJX = std::log(1.0 + (ATAU - 1.0) * DTAU600[L] / ATAU0X) / ATAULN;
            int JX = std::min(100, std::max(0, static_cast<int>(AJX + 0.5)));
            JXTRA[L] = JX;
            ATAU0X = ATAU0X * std::pow(ATAU, JX);
        }
    }

    // Check overflow of arrays: if total expanded levels * 2 > NX,
    // zero out JXTRA for this layer and all layers below.
    int JTOTL = L1X + 2;
    for (int L = L1X - 1; L >= 0; --L) {
        JTOTL += JXTRA[L];
        if (JTOTL * 2 > NX) {
            for (int LL = L; LL >= 0; --LL) {
                JXTRA[LL] = 0;
            }
            break;
        }
    }
}

// =========================================================================
// OPTICL — Liquid Water Cloud Optical Properties
// =========================================================================
/// Compute optical properties for liquid water clouds (cloud type 1) by
/// interpolating QCC/SCC/PCC tables against effective radius using the
/// RCC grid.
///
/// Ported from OPTICL in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param REFF   Effective radius of liquid water cloud (microns)
/// @param TEFF   Effective temperature (K) — unused for liquid clouds
/// @param DDENS  Output: density of cloud particle (g/cm³)
/// @param QQEXT  Output: extinction efficiency per S-bin [S_]
/// @param SSALB  Output: single-scattering albedo per S-bin [S_]
/// @param SSLEG  Output: phase function Legendre coefficients [8*S_],
///               accessed as SSLEG[L + 8*J] for L=0..7, J=0..S_-1
/// @param state  CloudJState containing cloud scattering tables
inline void OPTICL(double REFF, double TEFF, double& DDENS,
                   double* QQEXT, double* SSALB, double* SSLEG,
                   const CloudJState& state)
{
    (void)TEFF;  // unused for liquid water clouds

    constexpr int K = 0;  // liquid water (Fortran K=1 → C++ K=0)

    DDENS = state.DCC[K];

    // Find interpolation bracket in RCC grid
    // Fortran: I=1, loop NR=2..MCC-1, if REFF > RCC(NR,K) then I=NR
    // C++ 0-based: I=0, loop NR=1..MCC-2, if REFF > RCC[NR][K] then I=NR
    int I = 0;
    for (int NR = 1; NR <= state.MCC - 2; ++NR) {
        if (REFF > state.RCC[NR][K]) {
            I = NR;
        }
    }

    // Fractional position between RCC[I] and RCC[I+1]
    double FNR = (REFF - state.RCC[I][K]) / (state.RCC[I + 1][K] - state.RCC[I][K]);
    FNR = std::min(1.0, std::max(0.0, FNR));

    // Interpolate optical properties for each S-bin
    for (int J = 0; J < S_; ++J) {
        QQEXT[J] = state.QCC[J][I][K] + FNR * (state.QCC[J][I + 1][K] - state.QCC[J][I][K]);
        SSALB[J] = state.SCC[J][I][K] + FNR * (state.SCC[J][I + 1][K] - state.SCC[J][I][K]);
        for (int L = 0; L < 8; ++L) {
            SSLEG[L + 8 * J] = state.PCC[L][J][I][K] +
                               FNR * (state.PCC[L][J][I + 1][K] - state.PCC[L][J][I][K]);
        }
    }
}

// =========================================================================
// OPTICI — Ice Cloud Optical Properties
// =========================================================================
/// Compute optical properties for ice clouds by interpolating QCC/SCC/PCC
/// tables against effective radius. Selects between irregular-ice (cloud
/// type 2, warm: T >= 233.15K) and hexagonal-ice (cloud type 3, cold:
/// T < 233.15K).
///
/// Ported from OPTICI in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param REFF   Effective radius of ice cloud (microns)
/// @param TEFF   Effective temperature (K) — selects ice crystal type
/// @param DDENS  Output: density of cloud particle (g/cm³)
/// @param QQEXT  Output: extinction efficiency per S-bin [S_]
/// @param SSALB  Output: single-scattering albedo per S-bin [S_]
/// @param SSLEG  Output: phase function Legendre coefficients [8*S_],
///               accessed as SSLEG[L + 8*J] for L=0..7, J=0..S_-1
/// @param state  CloudJState containing cloud scattering tables
inline void OPTICI(double REFF, double TEFF, double& DDENS,
                   double* QQEXT, double* SSALB, double* SSLEG,
                   const CloudJState& state)
{
    // Select ice cloud type based on temperature
    // Fortran K=2 (irreg-ice, warm) → C++ K=1
    // Fortran K=3 (hexag-ice, cold) → C++ K=2
    int K;
    if (TEFF >= 233.15) {
        K = 1;  // ice irregular (warm)
    } else {
        K = 2;  // ice hexagonal (cold)
    }

    DDENS = state.DCC[K];

    // Find interpolation bracket in RCC grid
    int I = 0;
    for (int NR = 1; NR <= state.MCC - 2; ++NR) {
        if (REFF > state.RCC[NR][K]) {
            I = NR;
        }
    }

    // Fractional position between RCC[I] and RCC[I+1]
    double FNR = (REFF - state.RCC[I][K]) / (state.RCC[I + 1][K] - state.RCC[I][K]);
    FNR = std::min(1.0, std::max(0.0, FNR));

    // Interpolate optical properties for each S-bin
    for (int J = 0; J < S_; ++J) {
        QQEXT[J] = state.QCC[J][I][K] + FNR * (state.QCC[J][I + 1][K] - state.QCC[J][I][K]);
        SSALB[J] = state.SCC[J][I][K] + FNR * (state.SCC[J][I + 1][K] - state.SCC[J][I][K]);
        for (int L = 0; L < 8; ++L) {
            SSLEG[L + 8 * J] = state.PCC[L][J][I][K] +
                               FNR * (state.PCC[L][J][I + 1][K] - state.PCC[L][J][I][K]);
        }
    }
}

// =========================================================================
// OPTICS — UCI Stratospheric Sulfate Aerosol Optical Properties
// =========================================================================
/// Compute optical depth, single-scattering albedo, and phase function for
/// UCI stratospheric sulfate aerosol (SSA) using temperature-interpolated tables.
///
/// K=1 → background (KK=4), K=2 → volcanic (KK=13).
/// Extinction = 0.75 * Q / (Reff * density).
///
/// Ported from OPTICS in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param OPTD   Output optical depth per wavelength bin [S_]
/// @param SSALB  Output single-scattering albedo per bin [S_]
/// @param SLEG   Output phase function Legendre coefficients [8*S_], stored as SLEG[I + 8*J]
/// @param PATH   Path (g/m2) of aerosol
/// @param K      Input aerosol index (1 or 2); output: mapped KK index
/// @param state  CloudJState containing SSA tables (QSS, SSS, PSS, RSS, DSS)
/// @param rc     Return code (0=success, -1=error)
inline void OPTICS(double* OPTD, double* SSALB, double* SLEG,
                   double PATH, int& K, const CloudJState& state, int& rc)
{
    int KK;

    if (K == 1) {
        KK = 4 - 1;   // C++ 0-based: Fortran KK=4 → index 3
    } else if (K == 2) {
        KK = 13 - 1;  // C++ 0-based: Fortran KK=13 → index 12
    } else {
        CLOUDJ_ERROR("SSA index out-of-range",
                     "OPTICS in photo_jx.hpp", rc);
        return;
    }

    double REFF = state.RSS[KK];
    double RHO  = state.DSS[KK];

    for (int J = 0; J < S_; ++J) {
        // extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * density(g/cm3)]
        double XTINCT = 0.75 * state.QSS[J][KK] / (REFF * RHO);
        OPTD[J]  = PATH * XTINCT;
        SSALB[J] = state.SSS[J][KK];
        for (int I = 0; I < 8; ++I) {
            SLEG[I + 8 * J] = state.PSS[I][J][KK];
        }
    }

    K = KK + 1;  // Return Fortran 1-based KK to caller
}

// =========================================================================
// OPTICG — GeoMIP Stratospheric Sulfate Aerosol Optical Properties
// =========================================================================
/// Compute optical depth, SSA, and phase function for GeoMIP aerosol data.
/// K = 1001:1015 corresponds to different effective radii.
///
/// Ported from OPTICG in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param OPTD   Output optical depth per wavelength bin [S_]
/// @param SSALB  Output single-scattering albedo per bin [S_]
/// @param SLEG   Output phase function Legendre coefficients [8*S_], stored as SLEG[I + 8*J]
/// @param PATH   Path (g/m2) of aerosol
/// @param K      Aerosol index (1001..1000+NGG)
/// @param state  CloudJState containing GeoMIP tables (QGG, SGG, PGG, RGG, DGG)
inline void OPTICG(double* OPTD, double* SSALB, double* SLEG,
                   double PATH, int K, const CloudJState& state)
{
    // Fortran: KK = max(1, min(NGG, K-1000))
    // C++ 0-based: subtract 1 further
    int KK = std::max(1, std::min(state.NGG, K - 1000)) - 1;

    double REFF = state.RGG[KK];
    double RHO  = state.DGG[KK];

    for (int J = 0; J < S_; ++J) {
        // extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * density(g/cm3)]
        double XTINCT = 0.75 * state.QGG[J][KK] / (REFF * RHO);
        OPTD[J]  = PATH * XTINCT;
        SSALB[J] = state.SGG[J][KK];
        for (int I = 0; I < 8; ++I) {
            SLEG[I + 8 * J] = state.PGG[I][J][KK];
        }
    }
}

// =========================================================================
// OPTICA — Standard Aerosol Mie Optical Properties (Wavelength-Interpolated)
// =========================================================================
/// Compute optical depth, SSA, and phase function for standard aerosol types
/// using wavelength-interpolated Mie tables (5 reference wavelengths:
/// 200-300-400-600-999 nm). Interpolates using linear 1/wavelength mapping.
///
/// Ported from OPTICA in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param OPTD   Output optical depth per wavelength bin [S_]
/// @param SSALB  Output single-scattering albedo per bin [S_]
/// @param SLEG   Output phase function Legendre coefficients [8*S_], stored as SLEG[I + 8*J]
/// @param PATH   Path (g/m2) of aerosol
/// @param RELH   Relative humidity (0.00 to 1.00+) [unused in this routine, kept for interface]
/// @param K      Aerosol index (3..NAA); on output may be unchanged
/// @param state  CloudJState containing Mie tables (QAA, SAA, PAA, RAA, DAA, WL)
/// @param rc     Return code (0=success, -1=error)
inline void OPTICA(double* OPTD, double* SSALB, double* SLEG,
                   double PATH, double RELH, int& K,
                   const CloudJState& state, int& rc)
{
    (void)RELH;  // not used in this routine

    // K=1&2 are SSA values handled by OPTICS; K must be in [3..NAA]
    if (K > state.NAA || K < 3) {
        CLOUDJ_ERROR("Aerosol index out-of-range",
                     "OPTICA in photo_jx.hpp", rc);
        return;
    }

    // Convert to 0-based index for C++ arrays
    int KK = K - 1;

    double REFF = state.RAA[KK];
    double RHO  = state.DAA[KK];

    for (int J = 0; J < S_; ++J) {
        double WAVE = state.WL[J];  // mean wavelength (nm) for bin J
        int JMIE;
        double WAAX;
        double QAAX;

        // Pick pair of Mie wavelengths for interpolation
        JMIE = 0;  // 0-based index for wavelength pair starting at 200nm
        WAAX = (WAVE - 200.0) * 0.010;

        if (WAVE > 300.0) {
            JMIE = 1;
            WAAX = (WAVE - 300.0) * 0.010;
        }
        if (WAVE > 400.0) {
            JMIE = 2;
            WAAX = (WAVE - 400.0) * 0.005;
        }
        if (WAVE > 600.0) {
            JMIE = 3;
            WAAX = (WAVE - 600.0) * 0.0025;
        }

        if (WAVE > 999.0) {
            // Beyond last wavelength: Q scales as 1/wavelength
            QAAX = state.QAA[4][KK] * 999.0 / WAVE;
            SSALB[J] = state.SAA[4][KK];
            for (int I = 0; I < 8; ++I) {
                SLEG[I + 8 * J] = state.PAA[I][4][KK];
            }
        } else {
            // Interpolate between bracketing wavelengths
            WAAX = std::min(1.0, std::max(0.0, WAAX));
            QAAX = state.QAA[JMIE][KK] * (1.0 - WAAX) + state.QAA[JMIE + 1][KK] * WAAX;
            SSALB[J] = state.SAA[JMIE][KK] * (1.0 - WAAX) + state.SAA[JMIE + 1][KK] * WAAX;
            for (int I = 0; I < 8; ++I) {
                SLEG[I + 8 * J] = state.PAA[I][JMIE][KK] * (1.0 - WAAX)
                                + state.PAA[I][JMIE + 1][KK] * WAAX;
            }
        }

        // extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * density(g/cm3)]
        double XTINCT = 0.75 * QAAX / (REFF * RHO);
        OPTD[J] = PATH * XTINCT;
    }
}

// =========================================================================
// OPTICM — UMich Aerosol Optical Properties (Relative Humidity Interpolated)
// =========================================================================
/// Compute optical depth, SSA, and phase function for University of Michigan
/// aerosol types using the UMAER lookup table with relative humidity interpolation.
/// Phase function approximated as (2L+1)*g^L for L=0..7.
///
/// Ported from OPTICM in src/Core/cldj_fjx_sub_mod.F90.
///
/// UMAER dimensions: [3][6][21][33]
///   Index 0: SSA (single scattering albedo)
///   Index 1: g (asymmetry parameter)
///   Index 2: k-ext (extinction, m2/g)
///   J=0..5: wavelengths [200, 300, 400, (550,) 600, 1000 nm]
///   KR=0..20: relative humidity [0,5,10,...,90,95,99%]
///   L=0..32: UM aerosol types
///
/// @param OPTD   Output optical depth per wavelength bin [S_]
/// @param SSALB  Output single-scattering albedo per bin [S_]
/// @param SLEG   Output phase function Legendre coefficients [8*S_], stored as SLEG[I + 8*J]
/// @param PATH   Path (g/m2) of aerosol
/// @param RELH   Relative humidity (0.00 to 1.00)
/// @param K      Aerosol type index (1..33 in Fortran convention)
/// @param state  CloudJState containing UMAER table and WL
/// @param rc     Return code (0=success, -1=error)
inline void OPTICM(double* OPTD, double* SSALB, double* SLEG,
                   double PATH, double RELH, int K,
                   const CloudJState& state, int& rc)
{
    // Validate aerosol type index (Fortran 1-based: 1..33)
    int L = K;
    if (L < 1 || L > 33) {
        CLOUDJ_ERROR("Aerosol index out-of-range",
                     "OPTICM in photo_jx.hpp", rc);
        return;
    }

    // Convert to 0-based for C++ array access
    int LL = L - 1;

    // Pick nearest relative humidity bin (21 bins: 0..20)
    // Fortran: KR = 20*RELH + 1.5 (1-based), then clamp to [1,21]
    // C++ 0-based: KR = 20*RELH + 0.5, clamp to [0,20]
    int KR = static_cast<int>(20.0 * RELH + 0.5);
    KR = std::max(0, std::min(20, KR));

    for (int J = 0; J < S_; ++J) {
        double WAVE = state.WL[J];

        // Pick nearest Mie wavelength index (0-based)
        int JMIE = 0;  // use 200 nm for < 255 nm
        if (WAVE > 255.0) JMIE = 1;  // use 300 nm for 255-355 nm
        if (WAVE > 355.0) JMIE = 2;  // use 400 nm for 355-500 nm
        if (WAVE > 500.0) JMIE = 3;  // use 600 nm for 500-800 nm
        if (WAVE > 800.0) JMIE = 4;  // use 1000 nm for > 800 nm

        // UMAER[2][JMIE][KR][LL] = extinction (m2/g)
        double XTINCT = state.UMAER[2][JMIE][KR][LL];

        // Rescale/reduce optical depth as 1/wavelength for > 1000 nm
        if (WAVE > 1000.0) {
            XTINCT = XTINCT * 1000.0 / WAVE;
        }

        OPTD[J]  = PATH * XTINCT;
        SSALB[J] = state.UMAER[0][JMIE][KR][LL];

        // Asymmetry parameter g
        double GCOS = state.UMAER[1][JMIE][KR][LL];

        // Phase function: (2L+1) * g^L for L=0..7
        SLEG[0 + 8 * J] =  1.0;
        SLEG[1 + 8 * J] =  3.0 * GCOS;
        SLEG[2 + 8 * J] =  5.0 * GCOS * GCOS;
        SLEG[3 + 8 * J] =  7.0 * GCOS * GCOS * GCOS;
        SLEG[4 + 8 * J] =  9.0 * GCOS * GCOS * GCOS * GCOS;
        SLEG[5 + 8 * J] = 11.0 * GCOS * GCOS * GCOS * GCOS * GCOS;
        SLEG[6 + 8 * J] = 13.0 * GCOS * GCOS * GCOS * GCOS * GCOS * GCOS;
        SLEG[7 + 8 * J] = 15.0 * GCOS * GCOS * GCOS * GCOS * GCOS * GCOS * GCOS;
    }
}

// =========================================================================
// ACLIM_FJX — Load Fast-JX Climatology (T, O3, CH4) for Latitude & Month
// =========================================================================
/// Load fast-JX climatology – T & O3 & CH4 – for latitude & month on the
/// supplied pressure grid. Used to extend the atmospheric column above the
/// CTM top with reference profiles.
///
/// Ported from ACLIM_FJX in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param MONTH   Month index (1-12)
/// @param YLATD   Latitude in degrees (-90 to +90), NOT a pre-computed bin
///                index. The climatology bin is derived the same way as
///                Fortran: N = max(1, min(18, int(YLATD+99)/10)).
/// @param PPP     Edge-pressure of CTM layers (hPa), dimension [L1U+1]
/// @param TTT     Output: temperature profile (K), dimension [L1U]
/// @param OOO     Output: ozone mixing ratio (ppm), dimension [L1U]
/// @param CH4     Output: methane mixing ratio (ppb), dimension [L1U]
/// @param L1U     Number of CTM levels + 1 (dimension of CTM)
/// @param state   CloudJState containing T_REF, O_REF, H2O_REF, CH4_REF
inline void ACLIM_FJX(int MONTH, double YLATD,
    const double* PPP, double* TTT, double* OOO, double* CH4,
    int L1U, const CloudJState& state)
{
    double OREF2[LREF], TREF2[LREF], HREF2[LREF], CREF2[LREF];
    double PSTD[LREF + 1];

    // Select appropriate month (clamp to 1-12)
    int M = std::max(1, std::min(12, MONTH)) - 1;  // 0-based for C++ array

    // Select appropriate latitudinal profile.
    // Fortran: N = max(1, min(18, int(YLATD+99)/10)) using integer division.
    // Ported 1:1 (including Fortran's truncate-toward-zero int() then
    // integer division by 10) to match reference bin selection exactly.
    int N = std::max(1, std::min(JREF,
                static_cast<int>(YLATD + 99.0) / 10)) - 1;  // 0-based

    // Load reference profiles for selected lat/month
    for (int K = 0; K < LREF; ++K) {
        OREF2[K] = state.O_REF[K][N][M];
        TREF2[K] = state.T_REF[K][N][M];
        HREF2[K] = state.H2O_REF[K][N][M];
        CREF2[K] = state.CH4_REF[K][N][M];
    }

    // Set up pressure levels for O3/T climatology
    // Value given for each 2 km z* level applies from 1 km below to 1 km above
    PSTD[0] = std::max(PPP[0], 1000.0);
    PSTD[1] = 1000.0 * std::pow(10.0, -1.0 / 16.0);
    double DLOGP = std::pow(10.0, -2.0 / 16.0);
    for (int K = 2; K < LREF; ++K) {
        PSTD[K] = PSTD[K - 1] * DLOGP;
    }
    PSTD[LREF] = 0.0;

    // Apportion O3 and T onto CTM levels with mass (pressure) weighting
    for (int L = 0; L < L1U; ++L) {
        double F0 = 0.0;
        double T0 = 0.0;
        double H0 = 0.0;
        double C0 = 0.0;
        for (int K = 0; K < LREF; ++K) {
            double PC = std::min(PPP[L], PSTD[K]);
            double PB = std::max(PPP[L + 1], PSTD[K + 1]);
            if (PC > PB) {
                double XC = (PC - PB) / (PPP[L] - PPP[L + 1]);
                F0 += OREF2[K] * XC;
                T0 += TREF2[K] * XC;
                H0 += HREF2[K] * XC;
                C0 += CREF2[K] * XC;
            }
        }
        TTT[L] = T0;   // K
        OOO[L] = F0;   // ppm
        CH4[L] = C0;   // ppb
    }
}

// =========================================================================
// ACLIM_GEO — Compute Relative Humidity Profile
// =========================================================================
/// Calculates RH profile given mid-layer pressure (hPa), temperature (K),
/// and specific humidity (kg/kg). Uses Tetens formula for saturation vapor
/// pressure above freezing, and a modified form for ice below freezing.
///
/// Ported from ACLIM_RH in src/Core/cldj_fjx_sub_mod.F90.
/// (Note: renamed to ACLIM_GEO in C++ port per design specification)
///
/// @param PPP  Mid-layer pressure (hPa), dimension [L1U]
/// @param TTT  Temperature (K), dimension [L1U]
/// @param QQQ  Specific humidity (kg/kg), dimension [L1U]
/// @param RH   Output: relative humidity (0-1), dimension [L1U]
/// @param L1U  Number of CTM levels + 1
inline void ACLIM_GEO(const double* PPP, const double* TTT, const double* QQQ,
    double* RH, int L1U)
{
    const double eps = 287.04 / 461.50;

    for (int L = 0; L < L1U - 1; ++L) {
        double es;
        if (TTT[L] > 273.15) {
            // Tetens formula for liquid water
            double T = TTT[L] - 273.15;
            es = 6.112 * RadiativeSolver::exp_eval(17.67 * T / (T + 243.50));
        } else {
            // Ice saturation formula (Buck/Magnus form over ice)
            double T = TTT[L];  // in Kelvin
            double log_es = 23.33086 - 6111.72784 / T + 0.15215 * std::log(T);
            es = RadiativeSolver::exp_eval(log_es);
        }
        double qs = (eps * es) / (PPP[L] - es * (1.0 - eps));
        RH[L] = std::min(std::max(QQQ[L] / qs, 0.0), 1.0);
    }
    // Top layer: copy from layer below
    if (L1U > 1) {
        RH[L1U - 1] = RH[L1U - 2];
    }
}

// =========================================================================
// ACLIM_RH — Load GeoMIP SSA Climatology Profile
// =========================================================================
/// Load GEOMIP SSA climatology (vs P) for latitude & month given pressure grid.
/// Produces aerosol path (g/m2) and nearest effective radius index for each
/// model layer.
///
/// Ported from ACLIM_GEO in src/Core/cldj_fjx_sub_mod.F90.
/// (Note: renamed to ACLIM_RH in C++ port per design specification)
///
/// @param MONTH      Month index (1-12)
/// @param LATLOC     Latitude (degrees, -90 to +90)
/// @param AERS       Output: aerosol path per layer (g/m2), dimension [L1U]
/// @param NAER       Output: aerosol type index per layer (1001..1015), dimension [L1U]
/// @param PPP        Edge-pressure of CTM layers (hPa), dimension [L1U+1]
/// @param L1U        Number of CTM levels + 1
/// @param state      CloudJState containing GeoMIP data (P_GREF, R_GREF, X_GREF, RGG, NGG)
inline void ACLIM_RH(int MONTH, double YLATD, double* AERS, int* NAER,
    const double* PPP, int L1U, const CloudJState& state)
{
    double RREF2[LGREF + 2], XREF2[LGREF + 2], PREF2[LGREF + 2];
    double PSTD2[LGREF + 3];

    // Select appropriate month (clamp 1-12, convert to 0-based)
    int M = std::max(1, std::min(12, MONTH)) - 1;

    // Select appropriate latitudinal profile (1:64, delta = 2.7906)
    int N = 0;
    double YN = -86.5806;
    while (YLATD > YN) {
        N++;
        YN += 2.7906;
    }
    N = std::max(0, std::min(63, N));  // 0-based, clamp to [0,63]

    // Load GeoMIP reference profiles (reverse order from P_GREF)
    RREF2[0] = 0.0;
    XREF2[0] = 0.0;
    for (int K = 0; K < LGREF; ++K) {
        PREF2[K + 1] = state.P_GREF[LGREF - 1 - K];
        RREF2[K + 1] = state.R_GREF[N][LGREF - 1 - K][M];
        XREF2[K + 1] = state.X_GREF[N][LGREF - 1 - K][M];
    }
    PREF2[0] = PREF2[1] * 1.001;
    PREF2[LGREF + 1] = PREF2[LGREF] * 0.999;

    // Set PSTD2 to boundaries between the PREF2 points
    PSTD2[0] = std::max(PPP[0], PREF2[0]);
    for (int L = 1; L <= LGREF + 1; ++L) {
        PSTD2[L] = 0.5 * (PREF2[L - 1] + PREF2[L]);
    }
    PSTD2[LGREF + 2] = 0.0;

    // Compute REDGE boundaries for RGG bins
    double REDGE[GGA_];
    for (int I = 1; I < state.NGG - 1; ++I) {
        REDGE[I] = 0.5 * (state.RGG[I] + state.RGG[I - 1]);
    }

    // Integrate for pressure-weighted averages in each layer
    for (int L = 0; L < L1U; ++L) {
        double X0 = 0.0;
        double RX0 = 0.0;
        for (int K = 0; K < LGREF; ++K) {
            double PC = std::min(PPP[L], PSTD2[K]);
            double PB = std::max(PPP[L + 1], PSTD2[K + 1]);
            if (PC > PB) {
                double XC = (PC - PB) / (PPP[L] - PPP[L + 1]);
                X0 += XREF2[K] * XC;
                RX0 += RREF2[K] * XREF2[K] * XC;
            }
        }

        // Aerosol path (g/m2): rescale from microg-H2SO4/kg-air with 75 wt%
        AERS[L] = G100 * (PPP[L] - PPP[L + 1]) * X0 * 1.3333e-6;
        NAER[L] = 1001;

        // Pick nearest R-eff and rescale mass
        if (X0 > 0.0) {
            double R0 = RX0 / X0;
            int IGG = 0;
            for (int I = 1; I < state.NGG - 1; ++I) {
                if (R0 > REDGE[I]) {
                    IGG = I;
                }
            }
            NAER[L] = std::min(1000 + IGG + 1, 1000 + GGA_);  // +1 for 1-based
            AERS[L] = AERS[L] * state.RGG[IGG] / R0;
        }
    }
}

// =========================================================================
// JP_ATM0 — Print Formatted Atmospheric Column Diagnostic Summary
// =========================================================================
/// Print a short diagnostic summary of the atmospheric column.
/// Shows altitude, pressure, temperature, air density, ozone density,
/// and cumulative O2/O3 columns.
///
/// Ported from JP_ATM0 in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param PPJ   Edge-pressures (hPa), dimension [LU+2]
/// @param TTJ   Temperature (K) per layer, dimension [LU+1]
/// @param DDJ   Air column density per layer, dimension [LU+1]
/// @param OOJ   Ozone column density per layer, dimension [LU+1]
/// @param ZZJ   Altitude of layer edges (cm), dimension [LU+2]
/// @param LU    Number of CTM layers
inline void JP_ATM0(const double* PPJ, const double* TTJ, const double* DDJ,
    const double* OOJ, const double* ZZJ, int LU)
{
    std::printf("   L z(km)     p      T       d(air)   d(O3)"
                "  col(O2)  col(O3)      d-TAU   SS-alb"
                "  g(cos) CTM lyr=>\n");

    // Top edge (L = LU+2 in Fortran → index LU+1 in C++)
    int L_top = LU + 1;
    std::printf(" %3d%6.2f%10.3f\n",
        LU + 2, ZZJ[L_top] * 1.0e-5, PPJ[L_top]);

    double XCOLO2 = 0.0;
    double XCOLO3 = 0.0;
    double ZTOP = ZZJ[LU + 1];

    // Loop from LU+1 down to 1 (Fortran), i.e. index LU down to 0 (C++)
    for (int L = LU; L >= 0; --L) {
        XCOLO2 += DDJ[L] * 0.20948;
        XCOLO3 += OOJ[L];
        double DELZ = ZTOP - ZZJ[L];
        ZTOP = ZZJ[L];
        double ZKM = ZZJ[L] * 1.0e-5;
        double DAIR = DDJ[L] / DELZ;
        double DOZO = OOJ[L] / DELZ;
        std::printf(" %3d%6.2f%10.3f%7.2f%9.2e%9.2e%9.2e%9.2e\n",
            L + 1, ZKM, PPJ[L], TTJ[L], DAIR, DOZO, XCOLO2, XCOLO3);
    }
}

// =========================================================================
// JP_ATM — Print Per-Layer Optical Property Diagnostics
// =========================================================================
/// Print detailed per-layer diagnostics of the atmospheric column including
/// optical depth, single-scattering albedo, and asymmetry parameter.
/// Called when LPRTJ=true during PHOTO_JX execution.
///
/// Ported from JP_ATM in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param PPJ     Edge-pressures (hPa), dimension [LU+2]
/// @param TTJ     Temperature (K) per layer, dimension [LU+1]
/// @param DDJ     Air column density per layer, dimension [LU+1]
/// @param OOJ     Ozone column density per layer, dimension [LU+1]
/// @param HHJ     H2O column density per layer, dimension [LU+1]
/// @param ZZJ     Altitude of layer edges (cm), dimension [LU+2]
/// @param DTAU6   Optical depth per layer (600nm), dimension [LU+1]
/// @param POMEG6  Phase function moments per layer [8][LU+1], stored as POMEG6[I + 8*L]
/// @param JXTRA   Extra sub-layers per level, dimension [LU+1]
/// @param LU      Number of CTM layers
inline void JP_ATM(const double* PPJ, const double* TTJ, const double* DDJ,
    const double* OOJ, const double* HHJ, const double* ZZJ,
    const double* DTAU6, const double* POMEG6, const int* JXTRA, int LU)
{
    std::printf("   L z(km)     p      T   "
                "    d(air)   d(O3)  col(O2)  col(O3)  col(H2O)   d-TAU   SS-alb"
                "  g(cos) CTM lyr=>\n");

    // Top edge
    int L_top = LU + 1;
    std::printf(" %3d%6.2f%10.3f\n",
        LU + 2, ZZJ[L_top] * 1.0e-5, PPJ[L_top]);

    double XCOLO2 = 0.0;
    double XCOLO3 = 0.0;
    double XCOLH2O = 0.0;
    double ZTOP = ZZJ[LU + 1];

    // Loop from LU+1 down to 1 (Fortran), i.e. index LU down to 0 (C++)
    for (int L = LU; L >= 0; --L) {
        XCOLO2 += DDJ[L] * 0.20948;
        XCOLO3 += OOJ[L];
        XCOLH2O += HHJ[L];
        double DELZ = ZTOP - ZZJ[L];
        ZTOP = ZZJ[L];
        double ZKM = ZZJ[L] * 1.0e-5;
        double DAIR = DDJ[L] / DELZ;
        double DOZO = OOJ[L] / DELZ;
        // POMEG6 stored as [8*(LU+1)]: POMEG6(I,L) in Fortran = POMEG6[(I-1) + 8*(L-1)]
        // In C++ 0-based: POMEG6[0 + 8*L] = ssa, POMEG6[1 + 8*L]/3 = g
        double ssa = POMEG6[0 + 8 * L];
        double gcos = POMEG6[1 + 8 * L] / 3.0;
        std::printf(" %3d%6.2f%10.3f%7.2f%9.2e%9.2e%9.2e%9.2e%9.2e%10.4f%8.5f%8.5f%3d\n",
            L + 1, ZKM, PPJ[L], TTJ[L], DAIR, DOZO,
            XCOLO2, XCOLO3, XCOLH2O,
            DTAU6[L], ssa, gcos, JXTRA[L]);
    }
}

// =========================================================================
// PHOTO_JX scratch workspace — persistent, reused across calls
// =========================================================================
/// Holds all per-call scratch buffers PHOTO_JX previously heap-allocated
/// on every invocation. A single thread_local instance is reused across
/// calls so std::vector::resize() becomes a no-op in steady state (no
/// per-call allocation). thread_local (NOT static) is mandatory: with
/// CLOUDJ_USE_OPENMP the wavelength loop runs concurrently and CLOUD_JX
/// may call PHOTO_JX up to 4x per column — each thread must own its own
/// scratch, never share it.
///
/// Sizes are IDENTICAL to the previous local declarations. Buffers are
/// resized at the top of PHOTO_JX; resize() only reallocates when the
/// requested size differs from the current size.
struct PhotoJXWorkspace {
    // Column arrays (fully overwritten each call → no re-zero needed)
    std::vector<double> PPJ, ZZJ, TTJ, HHJ, DDJ, RRJ, OOJ, CCJ;
    // AMG: initialized to 1.0 then conditionally overwritten → must be
    // re-filled with 1.0 each call.
    std::vector<double> AMG;
    // AMF: zeroed internally by the SPHERE* routines → no external zero.
    std::vector<double> AMF_data;
    // Optical-depth accumulators (+=) → must be zeroed each call.
    std::vector<double> OD_arr, SSA_arr, SLEG_arr, OD600, FFXTAU;
    // DTAUX/POMEGAX: fully overwritten in the K/L transform loop before the
    // radiative solve, BUT only the S_ (=WW) active bins are written while
    // the buffer is sized L1U*WW; the remaining bins were 0 from the fresh
    // allocation. Zero to preserve that.
    std::vector<double> DTAUX_data, POMEGAX_data;
    std::vector<int> JXTRA; // fully written by EXTRAL1
    // OPMIE output buffers. flxd/flxd0 are accumulated (+=) inside OPMIE;
    // the rest are direct-assigned but zeroed here for safety.
    std::vector<double> AVGF_data, FJTOP_data, FJBOT_data, FIBOT_data,
        FSBOT_data, FJFLX_data, FLXD_data, FLXD0_data;
    std::vector<double> FFF_data; // fully overwritten for active bins; zeroed
    // Heating-rate accumulators (+=) → must be zeroed each call.
    std::vector<double> FFX, FFXNET;
    std::vector<double> FLXJ; // fully overwritten before read each K
    // JRATET flat output buffer (reused; JRATET zeroes what it needs)
    std::vector<double> VALJL_flat;
    // Persistent radiative-solver workspace (its internal buffers persist;
    // ws.resize(nd) is a cheap no-op when nd is unchanged).
    RadiativeSolver::Workspace ws;
    // SpecData rebuilt from CloudJState; state is const/unchanging so we
    // only rebuild when NJX changes (spec_njx == -1 means "not built yet").
    Photolysis::SpecData spec;
    int spec_njx = -1;
};

// =========================================================================
// PHOTO_JX — Main Column Photolysis Physics Driver
// =========================================================================
/// Gateway to single column fast-JX calculations.
/// Calculates J-values for a single Independent Column Atmosphere (ICA).
///
/// Ported from PHOTO_JX in src/Core/cldj_fjx_sub_mod.F90.
///
/// @param U0      cos(solar zenith angle)
/// @param SZA     solar zenith angle (degrees)
/// @param RFL     surface albedo [5][W_+W_r] (angles 1-4 + U0)
/// @param SOLF    solar flux scaling factor
/// @param LPRTJ   diagnostic print flag
/// @param PPP     pressure edges [L1U+1] (hPa)
/// @param ZZZ     altitude edges [L1U+1] (cm, geopotential)
/// @param TTT     temperature [L1U] (K)
/// @param HHH     H2O column density [L1U] (molecules/cm2)
/// @param DDD     dry-air column density [L1U] (molecules/cm2)
/// @param RRR     relative humidity [L1U]
/// @param OOO     O3 column density [L1U] (molecules/cm2)
/// @param CCC     CH4 column density [L1U]
/// @param LWP     liquid water path [L1U] (g/m2)
/// @param IWP     ice water path [L1U] (g/m2)
/// @param REFFL   liquid cloud effective radius [L1U] (microns)
/// @param REFFI   ice cloud effective radius [L1U] (microns)
/// @param AERSP   aerosol path [L1U*ANU] column-major (g/m2)
/// @param NDXAER  aerosol type index [L1U*ANU] column-major
/// @param L1U     number of levels (layers+1)
/// @param ANU     number of aerosol types per layer
/// @param NJXU    number of J-values requested
/// @param VALJXX  output J-values [(L1U-1)*NJXU] column-major
/// @param SKPERD  output heating rates [(S_+2)*L1U] column-major
/// @param SWMSQ   output energy budget [6]
/// @param OD18    output 600nm optical depth [L1U]
/// @param LDARK   output dark flag
/// @param FSBOT   output direct surface flux [W_+W_r]
/// @param FJBOT   output diffuse surface flux [W_+W_r]
/// @param FLXD    output solar flux deposited [L1U*(W_+W_r)] column-major
/// @param FJFLX   output diffuse flux [L1U*(W_+W_r)] column-major
/// @param state   CloudJState with all tables and configuration
/// @param rc      return code (0=success)
// PHOTO_JX_BODY_START
inline void PHOTO_JX(
    double U0, double SZA,
    const double* RFL_flat, double SOLF, bool LPRTJ,
    const double* PPP, const double* ZZZ,
    const double* TTT, const double* HHH,
    const double* DDD, const double* RRR,
    const double* OOO, const double* CCC,
    const double* LWP_in, const double* IWP_in,
    const double* REFFL, const double* REFFI,
    const double* AERSP, const int* NDXAER,
    int L1U, int ANU, int NJXU,
    double* VALJXX,
    double* SKPERD,
    double* SWMSQ,
    double* OD18,
    bool& LDARK,
    double* FSBOT, double* FJBOT,
    double* FLXD_out, double* FJFLX_out,
    const CloudJState& state,
    int& rc)
{
    // --- Initialize outputs ---
    const int LU = L1U - 1;
    const int WW = W_ + W_r;  // total wavelength bins (=18 for v8.0)

    std::memset(VALJXX, 0, sizeof(double) * LU * NJXU);
    std::memset(SKPERD, 0, sizeof(double) * (S_ + 2) * L1U);
    std::memset(SWMSQ, 0, sizeof(double) * 6);
    std::memset(OD18, 0, sizeof(double) * L1U);
    LDARK = false;

    // --- 1. Check for dark conditions: SZA > 98 deg ---
    if (SZA > 98.0) {
        LDARK = true;
        return;
    }

    // --- Persistent per-thread scratch workspace (allocation-free steady
    //     state). thread_local, never shared across concurrent threads. ---
    thread_local PhotoJXWorkspace wsx;

    // Resize all scratch buffers to the sizes this call needs. resize() is a
    // no-op (no allocation) when the size is unchanged. Sizes below are
    // IDENTICAL to the previous per-call local declarations.
    const int WW_ = W_ + W_r;              // = WW below
    const int AMF_dim_ = L1U + 1;
    wsx.PPJ.resize(L1U + 1);
    wsx.ZZJ.resize(L1U + 1);
    wsx.TTJ.resize(L1U);
    wsx.HHJ.resize(L1U);
    wsx.DDJ.resize(L1U);
    wsx.RRJ.resize(L1U);
    wsx.OOJ.resize(L1U);
    wsx.CCJ.resize(L1U);
    wsx.AMG.resize(L1U);
    wsx.AMF_data.resize(AMF_dim_ * AMF_dim_);
    wsx.OD_arr.resize(S_ * L1U);
    wsx.SSA_arr.resize(S_ * L1U);
    wsx.SLEG_arr.resize(8 * S_ * L1U);
    wsx.OD600.resize(L1U);
    wsx.FFXTAU.resize(S_ * 4);
    wsx.DTAUX_data.resize(L1U * WW_);
    wsx.POMEGAX_data.resize(8 * L1U * WW_);
    wsx.JXTRA.resize(L1U);
    wsx.AVGF_data.resize(L1U * WW_);
    wsx.FJTOP_data.resize(WW_);
    wsx.FJBOT_data.resize(WW_);
    wsx.FIBOT_data.resize(5 * WW_);
    wsx.FSBOT_data.resize(WW_);
    wsx.FJFLX_data.resize(L1U * WW_);
    wsx.FLXD_data.resize(L1U * WW_);
    wsx.FLXD0_data.resize(WW_);
    wsx.FFF_data.resize(W_ * L1U);
    wsx.FFX.resize(S_ * L1U);
    wsx.FFXNET.resize(S_ * 8);
    wsx.FLXJ.resize(L1U);
    wsx.VALJL_flat.resize(LU * NJXU);

    // --- 2. Load column arrays ---
    // PPJ/ZZJ/TTJ/HHJ/DDJ/RRJ/OOJ/CCJ are FULLY WRITTEN in the loop below
    // (all indices) before any read → no explicit zeroing required.
    std::vector<double>& PPJ = wsx.PPJ;
    std::vector<double>& ZZJ = wsx.ZZJ;
    std::vector<double>& TTJ = wsx.TTJ;
    std::vector<double>& HHJ = wsx.HHJ;
    std::vector<double>& DDJ = wsx.DDJ;
    std::vector<double>& RRJ = wsx.RRJ;
    std::vector<double>& OOJ = wsx.OOJ;
    std::vector<double>& CCJ = wsx.CCJ;

    for (int L = 0; L < L1U; ++L) {
        PPJ[L] = PPP[L];
        ZZJ[L] = ZZZ[L];
        TTJ[L] = TTT[L];
        HHJ[L] = HHH[L];
        DDJ[L] = DDD[L];
        RRJ[L] = RRR[L];
        OOJ[L] = OOO[L];
        CCJ[L] = CCC[L];
    }
    PPJ[L1U] = 0.0;
    ZZJ[L1U] = ZZZ[L1U - 1] + ZZHT;
    // PHOTO_JX_PART2_PLACEHOLDER

    // --- 3. Convert geopotential→geometric if ATM0 >= 3; compute AMG ---
    // AMG: previously `std::vector<double> AMG(L1U, 1.0)`. Initialized to 1.0
    // then only conditionally overwritten (ATM0 >= 3) → MUST be re-filled with
    // 1.0 each call to preserve behavior.
    std::vector<double>& AMG = wsx.AMG;
    std::fill(AMG.begin(), AMG.end(), 1.0);

    if (state.ATM0 >= 3) {
        // Convert geopotential to geometric heights
        for (int L = 1; L <= L1U; ++L) {
            ZZJ[L] = ZZJ[L] / (1.0 - ZZJ[L] / RAD);
        }
        // Scale factor for area expansion & mass increase
        for (int L = 0; L < L1U; ++L) {
            double ZMID = 0.5 * (ZZJ[L] + ZZJ[L + 1]);
            AMG[L] = (1.0 + ZMID / RAD) * (1.0 + ZMID / RAD);
        }
    }

    // --- 4. Compute Air Mass Factors (AMF) ---
    // AMF_data is zeroed internally by each SPHERE* routine (verified: SPHERE1F
    // loops over all dim*dim entries setting AMF[i]=0.0) → no external zero.
    const int AMF_dim = AMF_dim_;
    std::vector<double>& AMF_data = wsx.AMF_data;

    if (state.ATM0 == 0) {
        SPHERE1F(U0, RAD, ZZJ.data(), ZZHT, AMF_data.data(), L1U);
    } else if (state.ATM0 == 1) {
        SPHERE1N(U0, RAD, ZZJ.data(), ZZHT, AMF_data.data(), L1U);
    } else {
        SPHERE1R(U0, RAD, ZZJ.data(), ZZHT, AMF_data.data(), L1U);
    }
    // PHOTO_JX_PART3_PLACEHOLDER

    // --- 5. Per-layer optical depth accumulation ---
    // Local scattering arrays: OD[S_][L1U], SSA[S_][L1U], SLEG[8][S_][L1U]
    // ACCUMULATED INTO (+=) across cloud/aerosol/gas contributions → MUST be
    // zeroed each call. OD600 is set to 0 per-L then += → zero for safety.
    // FFXTAU accumulated (+=) → zero each call.
    std::vector<double>& OD_arr = wsx.OD_arr;
    std::vector<double>& SSA_arr = wsx.SSA_arr;
    std::vector<double>& SLEG_arr = wsx.SLEG_arr;
    std::vector<double>& OD600 = wsx.OD600;
    std::vector<double>& FFXTAU = wsx.FFXTAU;
    std::fill(OD_arr.begin(), OD_arr.end(), 0.0);
    std::fill(SSA_arr.begin(), SSA_arr.end(), 0.0);
    std::fill(SLEG_arr.begin(), SLEG_arr.end(), 0.0);
    std::fill(OD600.begin(), OD600.end(), 0.0);
    std::fill(FFXTAU.begin(), FFXTAU.end(), 0.0);

    // Access macros for column-major arrays
    // OD(K,L) -> OD_arr[K + S_*L], SSA(K,L) -> SSA_arr[K + S_*L]
    // SLEG(I,K,L) -> SLEG_arr[I + 8*(K + S_*L)]
    #define OD_IDX(K,L)     ((K) + S_*(L))
    #define SSA_IDX(K,L)    ((K) + S_*(L))
    #define SLEG_IDX(I,K,L) ((I) + 8*((K) + S_*(L)))
    #define FFXTAU_IDX(K,J) ((K) + S_*(J))

    // Temporary arrays for OPTIC* routines
    double QQEXT[S_], SSALB_tmp[S_], SSLEG_tmp[8 * S_];
    double OPTX[S_];
    double DDENS;

    for (int L = 0; L < L1U; ++L) {
        OD600[L] = 0.0;

        // Initialize with Rayleigh scattering
        for (int K = 0; K < S_; ++K) {
            double ODRAY = DDJ[L] * state.QRAYL[K];
            OD_arr[OD_IDX(K, L)] = ODRAY;
            SSA_arr[SSA_IDX(K, L)] = ODRAY;
            SLEG_arr[SLEG_IDX(0, K, L)] = 1.0 * ODRAY;
            SLEG_arr[SLEG_IDX(2, K, L)] = 0.5 * ODRAY;
        }
    // PHOTO_JX_PART4_PLACEHOLDER

        // --- Liquid Water Cloud ---
        if (LWP_in[L] > 1.0e-5 && REFFL[L] > 0.1) {
            double RE_LIQ = REFFL[L];
            double TE_ICE = TTT[L];
            OPTICL(RE_LIQ, TE_ICE, DDENS, QQEXT, SSALB_tmp, SSLEG_tmp, state);

            for (int K = 0; K < S_; ++K) {
                double ODL = LWP_in[L] * 0.75 * QQEXT[K] / (RE_LIQ * DDENS);
                OD_arr[OD_IDX(K, L)] += ODL;
                SSA_arr[SSA_IDX(K, L)] += SSALB_tmp[K] * ODL;
                FFXTAU[FFXTAU_IDX(K, 2)] += ODL * (1.0 - SSALB_tmp[K]);
                FFXTAU[FFXTAU_IDX(K, 3)] += ODL;
                for (int I = 0; I < 8; ++I) {
                    SLEG_arr[SLEG_IDX(I, K, L)] += SSLEG_tmp[I + 8*K] * SSALB_tmp[K] * ODL;
                }
                if (K == S_ - 1) {  // bin 18 (0-based: S_-1)
                    OD600[L] += ODL;
                }
            }
        }

        // --- Ice Water Cloud ---
        if (IWP_in[L] > 1.0e-5 && REFFI[L] > 0.1) {
            double RE_ICE = REFFI[L];
            double TE_ICE = TTT[L];
            OPTICI(RE_ICE, TE_ICE, DDENS, QQEXT, SSALB_tmp, SSLEG_tmp, state);

            for (int K = 0; K < S_; ++K) {
                double ODL = IWP_in[L] * 0.75 * QQEXT[K] / (RE_ICE * DDENS);
                OD_arr[OD_IDX(K, L)] += ODL;
                SSA_arr[SSA_IDX(K, L)] += SSALB_tmp[K] * ODL;
                FFXTAU[FFXTAU_IDX(K, 2)] += ODL * (1.0 - SSALB_tmp[K]);
                FFXTAU[FFXTAU_IDX(K, 3)] += ODL;
                for (int I = 0; I < 8; ++I) {
                    SLEG_arr[SLEG_IDX(I, K, L)] += SSLEG_tmp[I + 8*K] * SSALB_tmp[K] * ODL;
                }
                if (K == S_ - 1) {
                    OD600[L] += ODL;
                }
            }
        }
    // PHOTO_JX_PART5_PLACEHOLDER

        // --- Strat Sulfate Aerosol (SSA): index 1 or 2 ---
        for (int M = 0; M < ANU; ++M) {
            int NAER = NDXAER[L + L1U * M];  // column-major [L1U][ANU]
            if (NAER == 1 || NAER == 2) {
                double PATH = AERSP[L + L1U * M];
                if (PATH > 0.0) {
                    int NAER_copy = NAER;
                    OPTICS(OPTX, SSALB_tmp, SSLEG_tmp, PATH, NAER_copy, state, rc);
                    for (int K = 0; K < S_; ++K) {
                        OD_arr[OD_IDX(K, L)] += OPTX[K];
                        SSA_arr[SSA_IDX(K, L)] += SSALB_tmp[K] * OPTX[K];
                        for (int I = 0; I < 8; ++I) {
                            SLEG_arr[SLEG_IDX(I, K, L)] += SSLEG_tmp[I + 8*K] * SSALB_tmp[K] * OPTX[K];
                        }
                    }
                    OD600[L] += OPTX[S_ - 1];
                }
            }
        }

        // --- GeoMIP enhanced SSA: index > 1000 ---
        for (int M = 0; M < ANU; ++M) {
            int NAER = NDXAER[L + L1U * M];
            if (NAER > 1000) {
                double PATH = AERSP[L + L1U * M];
                if (PATH > 0.0) {
                    OPTICG(OPTX, SSALB_tmp, SSLEG_tmp, PATH, NAER, state);
                    for (int K = 0; K < S_; ++K) {
                        OD_arr[OD_IDX(K, L)] += OPTX[K];
                        SSA_arr[SSA_IDX(K, L)] += SSALB_tmp[K] * OPTX[K];
                        for (int I = 0; I < 8; ++I) {
                            SLEG_arr[SLEG_IDX(I, K, L)] += SSLEG_tmp[I + 8*K] * SSALB_tmp[K] * OPTX[K];
                        }
                    }
                    OD600[L] += OPTX[S_ - 1];
                }
            }
        }
    // PHOTO_JX_PART6_PLACEHOLDER

        // --- Standard aerosols: index 3..999 ---
        double RH = RRJ[L];
        for (int M = 0; M < ANU; ++M) {
            int NAER = NDXAER[L + L1U * M];
            double PATH = AERSP[L + L1U * M];
            if (PATH > 0.0 && NAER > 2 && NAER < 1000) {
                OPTICA(OPTX, SSALB_tmp, SSLEG_tmp, PATH, RH, NAER, state, rc);
                for (int K = 0; K < S_; ++K) {
                    OD_arr[OD_IDX(K, L)] += OPTX[K];
                    SSA_arr[SSA_IDX(K, L)] += SSALB_tmp[K] * OPTX[K];
                    for (int I = 0; I < 8; ++I) {
                        SLEG_arr[SLEG_IDX(I, K, L)] += SSLEG_tmp[I + 8*K] * SSALB_tmp[K] * OPTX[K];
                    }
                }
                OD600[L] += OPTX[S_ - 1];
            }
        }

        // --- UMich aerosols: index < 0 ---
        for (int M = 0; M < ANU; ++M) {
            int NAER = NDXAER[L + L1U * M];
            double PATH = AERSP[L + L1U * M];
            if (PATH > 0.0 && NAER < 0) {
                OPTICM(OPTX, SSALB_tmp, SSLEG_tmp, PATH, RH, -NAER, state, rc);
                for (int K = 0; K < S_; ++K) {
                    OD_arr[OD_IDX(K, L)] += OPTX[K];
                    SSA_arr[SSA_IDX(K, L)] += SSALB_tmp[K] * OPTX[K];
                    for (int I = 0; I < 8; ++I) {
                        SLEG_arr[SLEG_IDX(I, K, L)] += SSLEG_tmp[I + 8*K] * SSALB_tmp[K] * OPTX[K];
                    }
                }
                OD600[L] += OPTX[S_ - 1];
            }
        }
    // PHOTO_JX_PART7_PLACEHOLDER

        // --- Add H2O UV absorption (bins 1:W_) ---
        for (int K = 0; K < W_; ++K) {
            OD_arr[OD_IDX(K, L)] += state.QH2O[K] * HHJ[L];
        }

        // --- Add O2 & O3 absorption (bins 1:W_) ---
        for (int K = 0; K < W_; ++K) {
            double TTTX = TTJ[L];
            double XQO2 = CrossSections::interpolate(
                TTTX,
                state.TQQ[0][0], state.QO2[0][K],
                state.TQQ[1][0], state.QO2[1][K],
                state.TQQ[2][0], state.QO2[2][K],
                state.LQQ[0]);
            double XQO3 = CrossSections::interpolate(
                TTTX,
                state.TQQ[0][1], state.QO3[0][K],
                state.TQQ[1][1], state.QO3[1][K],
                state.TQQ[2][1], state.QO3[2][K],
                state.LQQ[1]);
            double ODABS = XQO3 * OOJ[L] + XQO2 * DDJ[L] * 0.20948;
            OD_arr[OD_IDX(K, L)] += ODABS;
        }

        // --- Renormalize SLEG by OD ---
        for (int K = 0; K < S_; ++K) {
            double od_kl = OD_arr[OD_IDX(K, L)];
            if (od_kl > 0.0) {
                for (int I = 0; I < 8; ++I) {
                    SLEG_arr[SLEG_IDX(I, K, L)] /= od_kl;
                }
            }
            FFXTAU[FFXTAU_IDX(K, 0)] += od_kl * (1.0 - SLEG_arr[SLEG_IDX(0, K, L)]);
            FFXTAU[FFXTAU_IDX(K, 1)] += od_kl;
        }

    } // end layer loop L
    // PHOTO_JX_PART8_PLACEHOLDER

    // --- 6. Transform OD/SLEG to DTAUX/POMEGAX format and copy OD600→OD18 ---
    // DTAUX(L,K) = OD(K,L), POMEGAX(I,L,K) = SLEG(I,K,L)
    // Only the S_ active bins are written (K in [0,S_)), while the buffers are
    // sized L1U*WW / 8*L1U*WW. WW == S_ for v8.0 so every element is written,
    // but zero for safety to preserve the fresh-allocation semantics.
    std::vector<double>& DTAUX_data = wsx.DTAUX_data;
    std::vector<double>& POMEGAX_data = wsx.POMEGAX_data;
    std::fill(DTAUX_data.begin(), DTAUX_data.end(), 0.0);
    std::fill(POMEGAX_data.begin(), POMEGAX_data.end(), 0.0);

    for (int K = 0; K < S_; ++K) {
        for (int L = 0; L < L1U; ++L) {
            DTAUX_data[L + L1U * K] = OD_arr[OD_IDX(K, L)];
            for (int I = 0; I < 8; ++I) {
                POMEGAX_data[I + 8 * (L + L1U * K)] = SLEG_arr[SLEG_IDX(I, K, L)];
            }
        }
    }

    for (int L = 0; L < L1U; ++L) {
        OD18[L] = OD600[L];
    }

    // --- 7. Call EXTRAL1 to determine sub-layer insertion ---
    // JXTRA is fully written by EXTRAL1 (all L1U entries) before use → no
    // explicit zeroing required. Zero anyway for safety/robustness.
    std::vector<int>& JXTRA = wsx.JXTRA;
    std::fill(JXTRA.begin(), JXTRA.end(), 0);
    EXTRAL1(OD600.data(), L1U, N_, state.ATAU, state.ATAU0, JXTRA.data());
    // PHOTO_JX_PART9_PLACEHOLDER

    // --- 8. Call OPMIE for actinic flux computation ---
    // Set up mdspan views for the radiative solver
    using RadiativeSolver::mdspan_2d;
    using RadiativeSolver::mdspan_2d_mut;
    using RadiativeSolver::mdspan_3d_mut;
    using RadiativeSolver::mdspan_1d;
    using RadiativeSolver::mdspan_1d_mut;

    mdspan_2d dtaux_view(DTAUX_data.data(), L1U, WW);
    mdspan_3d_mut pomegax_view(POMEGAX_data.data(), 8, L1U, WW);

    // AMF and AMG mdspan views
    mdspan_2d amf_view(AMF_data.data(), AMF_dim, AMF_dim);
    mdspan_1d amg_view(AMG.data(), L1U);

    // RFL view: the input is flat [5*(W_+W_r)], column-major
    mdspan_2d rfl_view(RFL_flat, 5, WW);

    // Output arrays for OPMIE.
    // OPMIE ACCUMULATES into flxd (flxd(l,k) += ...) and flxd0 (flxd0(k) +=
    // ...), so FLXD_data and FLXD0_data MUST be zeroed each call. The other
    // outputs (avgf/fjtop/fjbot/fibot/fsbot/fjflx) are direct-assigned by
    // OPMIE, but we zero all of them for safety (correctness first).
    std::vector<double>& AVGF_data = wsx.AVGF_data;
    std::vector<double>& FJTOP_data = wsx.FJTOP_data;
    std::vector<double>& FJBOT_data = wsx.FJBOT_data;
    std::vector<double>& FIBOT_data = wsx.FIBOT_data;
    std::vector<double>& FSBOT_data = wsx.FSBOT_data;
    std::vector<double>& FJFLX_data = wsx.FJFLX_data;
    std::vector<double>& FLXD_data = wsx.FLXD_data;
    std::vector<double>& FLXD0_data = wsx.FLXD0_data;
    std::fill(AVGF_data.begin(), AVGF_data.end(), 0.0);
    std::fill(FJTOP_data.begin(), FJTOP_data.end(), 0.0);
    std::fill(FJBOT_data.begin(), FJBOT_data.end(), 0.0);
    std::fill(FIBOT_data.begin(), FIBOT_data.end(), 0.0);
    std::fill(FSBOT_data.begin(), FSBOT_data.end(), 0.0);
    std::fill(FJFLX_data.begin(), FJFLX_data.end(), 0.0);
    std::fill(FLXD_data.begin(), FLXD_data.end(), 0.0);
    std::fill(FLXD0_data.begin(), FLXD0_data.end(), 0.0);

    mdspan_2d_mut avgf_view(AVGF_data.data(), L1U, WW);
    mdspan_1d_mut fjtop_view(FJTOP_data.data(), WW);
    mdspan_1d_mut fjbot_view(FJBOT_data.data(), WW);
    mdspan_2d_mut fibot_view(FIBOT_data.data(), 5, WW);
    mdspan_1d_mut fsbot_view(FSBOT_data.data(), WW);
    mdspan_2d_mut fjflx_view(FJFLX_data.data(), L1U, WW);
    mdspan_2d_mut flxd_view(FLXD_data.data(), L1U, WW);
    mdspan_1d_mut flxd0_view(FLXD0_data.data(), WW);

    // Workspace for OPMIE — hoisted into the thread_local scratch so its
    // internal buffers persist across calls; ws.resize(nd) is a cheap no-op
    // when nd is unchanged.
    RadiativeSolver::Workspace& ws = wsx.ws;
    int jaddto = 0;
    for (int l = 0; l < L1U; ++l) jaddto += JXTRA[l];
    int nd = 2 * L1U + 2 * jaddto + 1;
    ws.resize(nd);
    ws.resize_opmie(nd);

    RadiativeSolver::OPMIE(
        dtaux_view, pomegax_view, U0, rfl_view, amf_view, amg_view,
        JXTRA, avgf_view, fjtop_view, fjbot_view, fibot_view,
        fsbot_view, fjflx_view, flxd_view, flxd0_view, LU, state.LDOKR,
        state.ATAU, ws);
    // PHOTO_JX_PART10_PLACEHOLDER

    // --- 9. Compute FFF (actinic flux * solar * FL) and call JRATET ---

    // --- 9. Compute FFF (actinic flux * solar * FL) and call JRATET ---
    // FFF(K,L) = SOLF * FL(K) * AVGF(L,K) for active bins.
    // Only active bins (LDOKR>0) and L in [0,LU) are written; other elements
    // must remain 0 → zero each call.
    std::vector<double>& FFF_data = wsx.FFF_data;
    std::fill(FFF_data.begin(), FFF_data.end(), 0.0);
    double PREF1 = 0.0, PREF2 = 0.0;

    for (int K = 0; K < W_; ++K) {
        if (state.LDOKR[K] > 0) {
            for (int L = 0; L < LU; ++L) {
                FFF_data[K + W_ * L] = SOLF * state.FL[K] * avgf_view(L, K);
            }
            PREF1 += FSBOT_data[K] * SOLF * state.FL[K] * state.FPAR[K];
            PREF2 += FJBOT_data[K] * SOLF * state.FL[K] * state.FPAR[K];
        }
    }

    // Create mdspan for FFF for JRATET (layout_left: [W_][LU])
    // JRATET expects fff(k, l) with k=wavelength, l=layer
    mdspan_2d_mut fff_view(FFF_data.data(), W_, L1U);

    // Build SpecData from CloudJState for JRATET.
    // `state` is const and its cross-section tables do not change between
    // calls, so we build the flat SpecData once (or whenever NJX changes) and
    // reuse it from the thread_local workspace. Profiled at ~0.8%.
    Photolysis::SpecData& spec = wsx.spec;
    if (wsx.spec_njx != state.NJX) {
    spec.nw = W_;
    spec.ns = S_;
    spec.njx = state.NJX;
    spec.lqq.resize(state.NJX);
    spec.sqq.resize(state.NJX);
    spec.tqq.resize(state.NJX * 3);
    spec.inv_t12.resize(state.NJX, 0.0);
    spec.inv_t23.resize(state.NJX, 0.0);

    for (int j = 0; j < state.NJX; ++j) {
        spec.lqq[j] = state.LQQ[j];
        spec.sqq[j] = state.SQQ[j];
        spec.tqq[j * 3 + 0] = state.TQQ[0][j];
        spec.tqq[j * 3 + 1] = state.TQQ[1][j];
        spec.tqq[j * 3 + 2] = state.TQQ[2][j];
        double dt12 = state.TQQ[1][j] - state.TQQ[0][j];
        double dt23 = state.TQQ[2][j] - state.TQQ[1][j];
        spec.inv_t12[j] = (dt12 != 0.0) ? 1.0 / dt12 : 0.0;
        spec.inv_t23[j] = (dt23 != 0.0) ? 1.0 / dt23 : 0.0;
    }

    // Cross-sections, flat layout: qo2/qo3/q1d[k * 3 + t],
    // qqq[(k * 3 + t) * NJX + j]
    spec.qo2.resize(W_ * 3);
    spec.qo3.resize(W_ * 3);
    spec.q1d.resize(W_ * 3);
    spec.qqq.resize(W_ * 3 * state.NJX);
    for (int k = 0; k < W_; ++k) {
        spec.qo2[k * 3 + 0] = state.QO2[0][k];
        spec.qo2[k * 3 + 1] = state.QO2[1][k];
        spec.qo2[k * 3 + 2] = state.QO2[2][k];

        spec.qo3[k * 3 + 0] = state.QO3[0][k];
        spec.qo3[k * 3 + 1] = state.QO3[1][k];
        spec.qo3[k * 3 + 2] = state.QO3[2][k];

        spec.q1d[k * 3 + 0] = state.Q1D[0][k];
        spec.q1d[k * 3 + 1] = state.Q1D[1][k];
        spec.q1d[k * 3 + 2] = state.Q1D[2][k];

        for (int t = 0; t < 3; ++t) {
            for (int j = 0; j < state.NJX; ++j) {
                spec.qqq[(k * 3 + t) * state.NJX + j] = state.QQQ[k][t][j];
            }
        }
    }
        wsx.spec_njx = state.NJX;
    } // end SpecData (re)build guard
    // PHOTO_JX_PART11_PLACEHOLDER

    // Call JRATET to compute J-values into the reused flat buffer.
    // Flat layout valjl_flat[L*NJXU + J] (row-major, matches old [L][J]).
    std::vector<double>& VALJL = wsx.VALJL_flat;
    Photolysis::JRATET(PPJ, TTJ, fff_view, VALJL, spec, LU, NJXU, rc);
    if (rc != CLDJ_SUCCESS) return;

    // Copy VALJL to output VALJXX [LU][NJXU] column-major. This is the
    // IDENTICAL transpose as before: VALJXX[L + LU*J] <- valjl[L][J], now
    // sourced from the flat buffer valjl[L*NJXU + J] (same L, same J).
    for (int L = 0; L < LU; ++L) {
        for (int J = 0; J < NJXU; ++J) {
            VALJXX[L + LU * J] = VALJL[L * NJXU + J];
        }
    }

    // Diagnostic J-value table print (Fortran cldj_fjx_sub_mod.F90:781-784,
    // inside the LPRTJ block).  PHOTO_JX is called once per quick-column
    // approximation under CLDFLAG=7 with LPRTJ true only on the first call,
    // so this prints the first-QCA table, exactly as the reference does.
    if (LPRTJ) {
        std::printf(" Fast-J ----J-values----\n");
        std::printf(" L=  ");
        for (int J = 0; J < state.NJX; ++J) {
            // Fortran a6: pad/truncate the title to exactly 6 chars, then 3x.
            std::string title = state.TITLEJX[J].substr(0, 6);
            title.resize(6, ' ');
            std::printf("%s   ", title.c_str());
        }
        std::printf("\n");
        for (int L = LU - 1; L >= 0; --L) {
            std::printf("%3d", L + 1);
            for (int J = 0; J < state.NJX; ++J) {
                std::printf("%s", format_fortran_e9_2(VALJXX[L + LU * J]).c_str());
            }
            std::printf("\n");
        }
        std::fflush(stdout);
    }

    // --- 10. Compute heating rates and energy budget ---
    // FFX(K,L) accumulates fractional absorbed flux per super-bin per layer.
    // FFX and FFXNET are ACCUMULATED INTO (+=) → MUST be zeroed each call.
    // FLXJ is fully overwritten (FLXJ[0], FLXJ[1..LU-1], FLXJ[LU]) before it
    // is read within each KG iteration → no zeroing strictly needed, zeroed
    // for safety.
    std::vector<double>& FFX = wsx.FFX;
    std::vector<double>& FFXNET = wsx.FFXNET;
    std::vector<double>& FLXJ = wsx.FLXJ;
    std::fill(FFX.begin(), FFX.end(), 0.0);
    std::fill(FFXNET.begin(), FFXNET.end(), 0.0);
    std::fill(FLXJ.begin(), FLXJ.end(), 0.0);
    #define FFX_IDX(K,L)    ((K) + S_*(L))
    #define FFXNET_IDX(K,J) ((K) + S_*(J))

    double FREFI = 0.0, FREFL = 0.0, FREFS = 0.0;

    int KG = 0;
    for (int K = 0; K < S_; ++K) {
        for (int JG = 0; JG < state.NSJSUB[K]; ++JG) {
            if (KG >= WW) break;
            if (state.LDOKR[KG] > 0) {
                double SJSUB_KJG = state.SJSUB[K][JG];
                FREFI += FLXD0_data[KG] * SOLF * state.FW[K] * SJSUB_KJG;
                FREFL += FJTOP_data[KG] * SOLF * state.FW[K] * SJSUB_KJG;
                FREFS += SOLF * state.FW[K] * SJSUB_KJG;

                double FABOT = (FJBOT_data[KG] + FSBOT_data[KG]) - FIBOT_data[4 + 5*KG];
                double FXBOT = FSBOT_data[KG] - FABOT;

                // FLXJ: net diffuse flux divergence per layer
                FLXJ[0] = FJFLX_data[0 + L1U * KG] - FXBOT;
                for (int L = 1; L < LU; ++L) {
                    FLXJ[L] = FJFLX_data[L + L1U * KG] - FJFLX_data[(L-1) + L1U * KG];
                }
                FLXJ[LU] = FJTOP_data[KG] - FJFLX_data[(LU-1) + L1U * KG];

                double FFX0 = 0.0;
                for (int L = 0; L < L1U; ++L) {
                    double dep = (FLXD_data[L + L1U * KG] - FLXJ[L]) * SJSUB_KJG;
                    FFX0 += dep;
                    FFX[FFX_IDX(K, L)] += dep;
                }
    // PHOTO_JX_PART12_PLACEHOLDER

                // Accumulate FFXNET budget diagnostics
                FFXNET[FFXNET_IDX(K, 0)] += FLXD0_data[KG] * SJSUB_KJG;
                FFXNET[FFXNET_IDX(K, 1)] += FSBOT_data[KG] * SJSUB_KJG;
                FFXNET[FFXNET_IDX(K, 2)] += (FLXD0_data[KG] + FSBOT_data[KG]) * SJSUB_KJG;
                FFXNET[FFXNET_IDX(K, 3)] += FJTOP_data[KG] * SJSUB_KJG;
                FFXNET[FFXNET_IDX(K, 4)] += FFX0;
                FFXNET[FFXNET_IDX(K, 5)] += FABOT * SJSUB_KJG;
                FFXNET[FFXNET_IDX(K, 6)] += FSBOT_data[KG] * SJSUB_KJG;
                FFXNET[FFXNET_IDX(K, 7)] += FJBOT_data[KG] * SJSUB_KJG;
            }
            KG++;
        }
    }
    // PHOTO_JX_PART13_PLACEHOLDER

    // --- Compute heating rates K/day per super-bin per layer ---
    for (int L = 0; L < L1U; ++L) {
        double DPKL = HeatFac_ / (PPP[L] - PPP[L + 1]);
        for (int K = 0; K < S_; ++K) {
            // SKPERD[(S_+2)*L1U] column-major: SKPERD[K + (S_+2)*L]
            SKPERD[K + (S_ + 2) * L] = FFX[FFX_IDX(K, L)] * state.FW[K] * SOLF * DPKL;
        }
        // Sum UV (bins 0..W_-1) and vis/IR (bins W_..S_-1)
        double sum_uv = 0.0, sum_ir = 0.0;
        for (int K = 0; K < W_; ++K) {
            sum_uv += SKPERD[K + (S_ + 2) * L];
        }
        for (int K = W_; K < S_; ++K) {
            sum_ir += SKPERD[K + (S_ + 2) * L];
        }
        SKPERD[S_ + (S_ + 2) * L] = sum_uv;
        SKPERD[S_ + 1 + (S_ + 2) * L] = sum_ir;
    }

    // --- Compute SWMSQ energy budget [6] ---
    for (int K = 0; K < S_; ++K) {
        SWMSQ[0] += FFXNET[FFXNET_IDX(K, 2)] * SOLF * state.FW[K];
        SWMSQ[1] += FFXNET[FFXNET_IDX(K, 3)] * SOLF * state.FW[K];
        SWMSQ[2] += FFXNET[FFXNET_IDX(K, 4)] * SOLF * state.FW[K];
        SWMSQ[3] += FFXNET[FFXNET_IDX(K, 5)] * SOLF * state.FW[K];
    }
    for (int K = 0; K < W_; ++K) {
        SWMSQ[4] += FFXNET[FFXNET_IDX(K, 6)] * state.FL[K] * state.FPAR[K] * SOLF;
        SWMSQ[5] += FFXNET[FFXNET_IDX(K, 7)] * state.FL[K] * state.FPAR[K] * SOLF;
    }
    // PHOTO_JX_PART14_PLACEHOLDER

    // --- Copy surface/layer flux outputs ---
    for (int K = 0; K < WW; ++K) {
        FSBOT[K] = FSBOT_data[K];
        FJBOT[K] = FJBOT_data[K];
    }
    for (int L = 0; L < L1U; ++L) {
        for (int K = 0; K < WW; ++K) {
            FLXD_out[L + L1U * K] = FLXD_data[L + L1U * K];
            FJFLX_out[L + L1U * K] = FJFLX_data[L + L1U * K];
        }
    }

    // Clean up macros
    #undef OD_IDX
    #undef SSA_IDX
    #undef SLEG_IDX
    #undef FFXTAU_IDX
    #undef FFX_IDX
    #undef FFXNET_IDX

    (void)LPRTJ;  // Diagnostic printing not implemented in C++ port
}

} // namespace PhotoJX
} // namespace CloudJ

#endif // CLOUDJ_PHOTO_JX_HPP
