#ifndef CLOUDJ_CLOUD_JX_HPP
#define CLOUDJ_CLOUD_JX_HPP

#include <cloudj/state.hpp>
#include <cloudj/photo_jx.hpp>
#include <cloudj/error.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace CloudJ {
namespace CloudOverlap {

// =========================================================================
// HEAPSORT_A - Classic heapsort for sorting ICAs by optical depth
// =========================================================================
/// Sorts real array A[N] into ASCENDING order, placing sorted values in AX[N].
/// Returns indexing IX[N] that records the permutation:
///   A[IX[J]] ==> AX[J], so IX[0] = original location of smallest A
///   and IX[N-1] = original location of largest A.
///
/// Ported from HEAPSORT_A in src/Core/cldj_sub_mod.F90.
///
/// @param N   Number of elements to sort
/// @param A   Input array of values (dimension >= N), 0-indexed
/// @param AX  Output sorted array (dimension >= N), 0-indexed
/// @param IX  Output index permutation array (dimension >= N), 0-indexed
/// @param ND  Array dimension (>= N)
inline void HEAPSORT_A(int N, const double* A, double* AX, int* IX, int ND)
{
    (void)ND;  // used only for bounds in Fortran

    // Initialize: copy A into AX, set IX to identity
    for (int i = 0; i < N; ++i) {
        IX[i] = i;
        AX[i] = A[i];
    }

    if (N <= 1) return;

    int L  = N / 2;  // Fortran: L = N/2 + 1, but we pre-decrement below
    int IR = N - 1;  // last valid index (0-based)

    double RA;
    int IA;

    for (;;) {
        if (L > 0) {
            --L;
            RA = AX[L];
            IA = IX[L];
        } else {
            RA = AX[IR];
            IA = IX[IR];
            AX[IR] = AX[0];
            IX[IR] = IX[0];
            --IR;
            if (IR == 0) {
                AX[0] = RA;
                IX[0] = IA;
                return;
            }
        }
        int i = L;
        int j = L + L + 1;  // 0-based child: 2*i+1
        while (j <= IR) {
            if (j < IR) {
                if (AX[j] < AX[j + 1]) {
                    ++j;
                }
            }
            if (RA < AX[j]) {
                AX[i] = AX[j];
                IX[i] = IX[j];
                i = j;
                j = 2 * i + 1;  // 0-based child
            } else {
                break;
            }
        }
        AX[i] = RA;
        IX[i] = IA;
    }
}

// =========================================================================
// ICA_NR - Set up max-random overlap groups and quantize cloud fractions
// =========================================================================
/// Reads cloud fraction (CLDF), cloud OD (CLTAU), cloud index (CLDIW).
/// Derives max-correlated cloud overlaps. Sets up the ICAs.
///
/// NCLDF[L] has value 0:CBIN_ = quantized cloud fraction.
/// CFBIN[J] = cloud fraction assumed for bin J=0:CBIN_-1.
/// CLTAU[L] is readjusted for quantum bins to preserve CLDF*CLTAU.
///
/// Ported from ICA_NR in src/Core/cldj_sub_mod.F90.
///
/// @param CLDF    Cloud fraction per layer [LTOP], 0-indexed
/// @param CLTAU   Cloud optical depth per layer [LTOP], modified in-place
/// @param IWPX    Ice water path per layer [LTOP], modified in-place
/// @param LWPX    Liquid water path per layer [LTOP], modified in-place
/// @param ZZZ     Height of layer edges [LTOP+1], 0-indexed
/// @param CLDIW   Cloud type index per layer [LTOP], 0-indexed
/// @param LTOP    Number of cloud layers
/// @param CLDCOR  Cloud decorrelation parameter
/// @param LNRG    Number of max-overlap groups config (0, 3, or 6)
/// @param NCLDF   Output: quantized cloud fractions [LTOP], 0-indexed
/// @param CFBIN   Output: cloud fraction bin values [CBIN_]
/// @param GFNR    Output: cloud fraction quantum nos per group [9][CBIN_+1]
/// @param GCMX    Output: max cloud fraction bin per group [9]
/// @param GNR     Output: number of unique fractions per group [9]
/// @param GBOT    Output: bottom layer of each group [9]
/// @param GTOP    Output: top layer of each group [9]
/// @param GLVL    Output: group level index [9]
/// @param NRG     Output: number of overlap groups
/// @param NICA    Output: total number of ICAs
inline void ICA_NR(const double* CLDF, double* CLTAU, double* IWPX,
                   double* LWPX, const double* ZZZ, const int* CLDIW,
                   int LTOP, double CLDCOR, int LNRG,
                   int* NCLDF, double* CFBIN,
                   int GFNR[9][CBIN_ + 1], int* GCMX, int* GNR,
                   int* GBOT, int* GTOP, int* GLVL,
                   int& NRG, int& NICA)
{
    (void)CLDCOR;  // decorrelation applied in ICA_ALL, not here

    // Local parameters
    constexpr int NG_BRK = 0;
    constexpr int NRG6_  = 6;
    constexpr double Zbin[NRG6_] = {0.0e5, 1.5e5, 3.5e5, 6.0e5, 9.0e5, 13.0e5};

    // Initialize outputs
    NRG  = 0;
    NICA = 0;
    std::memset(NCLDF, 0, LTOP * sizeof(int));
    std::memset(GBOT,  0, 9 * sizeof(int));
    std::memset(GTOP,  0, 9 * sizeof(int));
    std::memset(GLVL,  0, 9 * sizeof(int));
    std::memset(GNR,   0, 9 * sizeof(int));
    std::memset(GCMX,  0, 9 * sizeof(int));
    for (int g = 0; g < 9; ++g)
        for (int c = 0; c <= CBIN_; ++c)
            GFNR[g][c] = 0;
    std::memset(CFBIN, 0, CBIN_ * sizeof(double));

    // Quantize cloud fractions into bins
    double FBIN = static_cast<double>(CBIN_);
    for (int K = 0; K < CBIN_; ++K) {
        CFBIN[K] = static_cast<double>(K + 1) / FBIN;
    }

    // Quantize cloud fractions and adjust to conserve TAU*CF
    constexpr double CLF_MIN = 0.02;
    constexpr double CLF_MAX = 1.0 - CLF_MIN;

    for (int L = 0; L < LTOP; ++L) {
        if (CLDF[L] < CLF_MIN || CLDIW[L] == 0) {
            NCLDF[L] = 0;
            CLTAU[L] = 0.0;
            IWPX[L]  = 0.0;
            LWPX[L]  = 0.0;
        } else if (CLDF[L] > CLF_MAX) {
            NCLDF[L] = static_cast<int>(FBIN);
        } else {
            double FSCALE2 = CLDF[L] * FBIN + 0.4999;
            NCLDF[L] = std::max(static_cast<int>(FSCALE2), 1);
        }
    }

    // Find cloud-top layer
    int LCLTOP = 0;
    for (int L = 0; L < LTOP; ++L) {
        if (NCLDF[L] > 0) {
            LCLTOP = L;
        }
    }
    if (LCLTOP == 0 && (LTOP == 0 || NCLDF[0] == 0)) {
        // No clouds at all
        NRG  = 0;
        NICA = 1;
        return;
    }
    // Make LCLTOP inclusive (0-based index of top cloud layer)
    // If only layer 0 has clouds, LCLTOP=0 is correct

    // Rescale LWPX, IWPX, CLTAU to preserve CLDF*TAU with quantized fractions
    for (int L = 0; L <= LCLTOP; ++L) {
        if (NCLDF[L] > 0) {
            double FSCALE = CLDF[L] / CFBIN[NCLDF[L] - 1];
            CLTAU[L] *= FSCALE;
            IWPX[L]  *= FSCALE;
            LWPX[L]  *= FSCALE;
        }
    }

    // Local group arrays
    int GMIN[9] = {0}, GMAX[9] = {0};

    if (LNRG == 0) {
        // ---- LNRG=0: Break by minimum cloud fraction threshold ----
        int L = 0;
        NRG = 0;
        while (L < LTOP) {
            if (NCLDF[L] > NG_BRK) {
                if (NRG >= 9) break;
                GMIN[NRG] = L;
                GMAX[NRG] = LCLTOP;
                for (int LL = L + 1; LL <= LCLTOP; ++LL) {
                    if (NCLDF[LL] <= NG_BRK) {
                        GMAX[NRG] = LL;
                        break;
                    }
                }
                L = GMAX[NRG] + 1;
                ++NRG;
            } else {
                ++L;
            }
        }
        NRG = std::max(NRG, 1);
        GMIN[0] = 0;
        GMAX[NRG - 1] = LTOP - 1;
        for (int N = 1; N < NRG; ++N) {
            GMIN[N] = GMAX[N - 1] + 1;
        }
        for (int N = 0; N < NRG; ++N) {
            GLVL[N] = N + 1;
        }

    } else if (LNRG == 3) {
        // ---- LNRG=3: Fixed 3 random-overlap groups ----
        int L1 = 0;
        int L2 = 8;  // Fortran L2=9, 0-based = 8
        int L3 = L2;
        for (int L = LCLTOP; L >= L2; --L) {
            if (CLDIW[L] == 1 || CLDIW[L] == 3) {
                L3 = L + 1;
                break;
            }
        }
        bool L1GRP = false, L2GRP = false, L3GRP = false;
        for (int L = L1; L < L2; ++L) {
            L1GRP = L1GRP || (NCLDF[L] > 0);
        }
        for (int L = L2; L < L3; ++L) {
            L2GRP = L2GRP || (NCLDF[L] > 0);
        }
        for (int L = L3; L <= LCLTOP; ++L) {
            L3GRP = L3GRP || (NCLDF[L] > 0);
        }
        NRG = 0;
        if (L1GRP) { GMIN[NRG] = L1; GMAX[NRG] = L2 - 1; ++NRG; }
        if (L2GRP) { GMIN[NRG] = L2; GMAX[NRG] = L3 - 1; ++NRG; }
        if (L3GRP) { GMIN[NRG] = L3; GMAX[NRG] = LCLTOP; ++NRG; }
        NRG = std::max(NRG, 1);
        GMIN[0] = 0;
        GMAX[NRG - 1] = LTOP - 1;
        for (int N = 0; N < NRG; ++N) {
            GLVL[N] = N + 1;
        }

    } else {
        // ---- LNRG=6 (default): Altitude-defined correlation groups ----
        // Find levels in each of the NRG6_ altitude-defined groups
        for (int N = 0; N < NRG6_ - 1; ++N) {
            GMAX[N] = 0;
        }
        for (int L = 0; L <= LCLTOP; ++L) {
            for (int N = 1; N < NRG6_; ++N) {
                if ((ZZZ[L] - ZZZ[0]) < Zbin[N]) {
                    GMAX[N - 1] = L;
                }
            }
        }
        GMIN[0] = 0;
        for (int N = 1; N < NRG6_; ++N) {
            GMIN[N] = GMAX[N - 1] + 1;
        }
        GMAX[NRG6_ - 1] = LCLTOP;

        // Find groups that contain clouds
        NRG = 0;
        int GMIN_tmp[9], GMAX_tmp[9];
        for (int N = 0; N < NRG6_; ++N) {
            bool L6GRP = false;
            for (int L = GMIN[N]; L <= GMAX[N]; ++L) {
                if (NCLDF[L] > 0) L6GRP = true;
            }
            if (L6GRP) {
                GMIN_tmp[NRG] = GMIN[N];
                GMAX_tmp[NRG] = GMAX[N];
                GLVL[NRG] = N + 1;  // 1-based group level
                ++NRG;
            }
        }
        for (int N = 0; N < NRG; ++N) {
            GMIN[N] = GMIN_tmp[N];
            GMAX[N] = GMAX_tmp[N];
        }

        // Pull off cirrus shields from top MAX-GRP as separate group
        if (NRG > 0) {
            int LCIRRUS = 0;
            for (int L = GMAX[NRG - 1]; L >= GMIN[NRG - 1]; --L) {
                if (NCLDF[L] > static_cast<int>(FBIN) / 2 && CLDIW[L] == 2) {
                    LCIRRUS = L;
                }
            }
            if (LCIRRUS > GMIN[NRG - 1]) {
                GMIN[NRG] = LCIRRUS;
                GMAX[NRG] = GMAX[NRG - 1];
                GLVL[NRG] = 7;
                GMAX[NRG - 1] = LCIRRUS - 1;
                ++NRG;
            }
        }
        // Avoid gaps
        NRG = std::max(NRG, 1);
        GMIN[0] = 0;
        GMAX[NRG - 1] = LTOP - 1;
        for (int N = 1; N < NRG; ++N) {
            GMIN[N] = GMAX[N - 1] + 1;
        }
    }

    // --- Finished selection of max-overlap groups ---
    // Assign GBOT/GTOP from GMIN/GMAX
    if (NRG == 0) {
        NRG = 1;
        GBOT[0] = 0;
        GTOP[0] = LTOP - 1;
    } else {
        GBOT[0] = 0;
        GTOP[0] = GMAX[0];
        for (int N = 1; N < NRG; ++N) {
            GBOT[N] = std::max(GTOP[N - 1] + 1, GMIN[N]);
            GTOP[N] = GMAX[N];
        }
        GTOP[NRG - 1] = LTOP - 1;
    }

    // For each max-overlap group, calculate unique cloud fractions
    int NSAME[CBIN_];
    for (int N = 0; N < NRG; ++N) {
        std::memset(NSAME, 0, sizeof(NSAME));
        GCMX[N] = 0;
        for (int L = GBOT[N]; L <= GTOP[N]; ++L) {
            if (NCLDF[L] > 0) {
                NSAME[NCLDF[L] - 1] = 1;  // 0-based bin index
                GCMX[N] = std::max(GCMX[N], NCLDF[L]);
            }
        }
        // Sort cloud fractions in decreasing order
        // Largest bin = CBIN_ (treated as 100%)
        GFNR[N][0] = CBIN_;
        int NC = 1;
        for (int I = CBIN_ - 2; I >= 0; --I) {
            if (NSAME[I] > 0) {
                GFNR[N][NC] = I + 1;  // 1-based bin value
                ++NC;
            }
        }
        GNR[N] = NC;
        GFNR[N][NC] = 0;
    }

    // Calculate total number of ICAs
    NICA = 1;
    int NICAX = 1;
    int NRGX = 0;
    for (int N = 0; N < NRG; ++N) {
        NICA = NICA * GNR[N];
        if (NICA <= ICA_) {
            NICAX = NICA;
            NRGX = N + 1;
        }
    }
    if (NICA > ICA_) {
        NICA = NICAX;
        NRG = NRGX;
    }
}

// =========================================================================
// ICA_ALL - Generate weights and total OD for all ICAs
// =========================================================================
/// Using the max-ran cloud overlap info from ICA_NR, generates all ICAs.
///   OCOL[i] = cloud optical depth (total) in ICA i
///   WCOL[i] = weight (fractional area) of ICA i
///
/// Ported from ICA_ALL in src/Core/cldj_sub_mod.F90.
///
/// @param CLF     Cloud fraction per layer [LTOP], 0-indexed
/// @param CLT     Cloud optical depth per layer [LTOP], 0-indexed
/// @param LTOP    Number of cloud layers
/// @param CFBIN   Cloud fraction bin values [CBIN_]
/// @param CLDCOR  Cloud decorrelation parameter
/// @param NCLDF   Quantized cloud fractions [LTOP]
/// @param GFNR    Cloud fraction quantum nos per group [9][CBIN_+1]
/// @param GCMX    Max cloud fraction bin per group [9]
/// @param GNR     Number of unique fractions per group [9]
/// @param GBOT    Bottom layer of each group [9]
/// @param GTOP    Top layer of each group [9]
/// @param GLVL    Group level index [9]
/// @param NRG     Number of overlap groups
/// @param NICA    Total number of ICAs
/// @param WCOL    Output: weight of each ICA [ICA_]
/// @param OCOL    Output: total OD of each ICA [ICA_]
inline void ICA_ALL(const double* CLF, const double* CLT, int LTOP,
                    const double* CFBIN, double CLDCOR,
                    const int* NCLDF,
                    const int GFNR[9][CBIN_ + 1], const int* GCMX,
                    const int* GNR, const int* GBOT, const int* GTOP,
                    const int* GLVL, int NRG, int NICA,
                    double* WCOL, double* OCOL)
{
    (void)CLF;   // cloud fractions already encoded in NCLDF
    (void)LTOP;  // layer range given by GBOT/GTOP

    // Initialize outputs
    for (int i = 0; i < NICA; ++i) {
        WCOL[i] = 0.0;
        OCOL[i] = 0.0;
    }

    // CF0: cloud fraction boundaries (0-based: CF0[0]=0, CF0[1..CBIN_]=CFBIN)
    double CF0[CBIN_ + 2];
    CF0[0] = 0.0;
    for (int L = 0; L < CBIN_; ++L) {
        CF0[L + 1] = CFBIN[L];
    }

    // Per-group data
    double FCMX[10];
    bool LGR_CLR[10];
    int GCLDY[10];
    double FWT[10][CBIN_ + 1];
    double FWTC[10][CBIN_ + 1];
    double FWTCC[10][CBIN_ + 1];

    for (int G = 0; G < NRG; ++G) {
        FCMX[G] = CF0[GCMX[G]];  // max cloud-fraction in MAX-GRP
        if (FCMX[G] < 0.99) {
            LGR_CLR[G] = true;   // 1st member of MAX-GRP G = clear sky
            GCLDY[G] = 2;
        } else {
            LGR_CLR[G] = false;
            GCLDY[G] = 1;
        }
        for (int I = 0; I < GNR[G]; ++I) {
            FWT[G][I]  = CF0[GFNR[G][I]] - CF0[GFNR[G][I + 1]];
            FWTC[G][I] = FWT[G][I];
            FWTCC[G][I] = FWT[G][I];
        }
    }
    if (NRG < 10) FCMX[NRG] = 0.0;

    // Pre-calculate correlation factors between adjacent groups
    for (int G = 0; G < NRG - 1; ++G) {
        bool LSKIP = (GCMX[G + 1] == 0 || GCMX[G + 1] == CBIN_
                   || GCMX[G] == 0 || GCMX[G] == CBIN_);
        if (!LSKIP) {
            double FIG2 = FCMX[G + 1];     // cloudy fract of upper group
            int GRP2    = GLVL[G + 1];      // upper G6 group number
            double FIG1 = FCMX[G];          // cloudy fract of current group
            int GRP1    = GLVL[G];          // current G6 group number
            double CORRFAC = std::pow(CLDCOR, GRP2 - GRP1);

            double GCORR = std::min({
                1.0 + CORRFAC * (1.0 / FIG2 - 1.0),
                1.0 / FIG2,
                1.0 / FIG1
            });

            for (int I = 1; I < GNR[G]; ++I) {
                FWTC[G][I]  = GCORR * FWT[G][I];
                FWTCC[G][I] = FWT[G][I] * (1.0 - GCORR * FIG2) / (1.0 - FIG2);
            }
            FWTC[G][0]  = 1.0 - FIG1 * GCORR;
            FWTCC[G][0] = 1.0 - FIG1 * (1.0 - GCORR * FIG2) / (1.0 - FIG2);
        }
    }

    // Generate all ICAs
    int IGNR[10];
    for (int I = 0; I < NICA; ++I) {
        double WTCOL = 1.0;
        double ODCOL = 0.0;

        // For each ICA, locate members of each group
        int II = I + 1;  // 1-based ICA index for modular arithmetic
        for (int G = 0; G < NRG; ++G) {
            IGNR[G] = ((II - 1) % GNR[G]);  // 0-based member index
            II = (II - 1) / GNR[G] + 1;
        }

        for (int G = 0; G < NRG; ++G) {
            int IG1 = IGNR[G];
            bool L_CLR1 = (GFNR[G][IG1] > GCMX[G]);

            bool L_CLR2;
            if (G == NRG - 1) {
                L_CLR2 = true;
            } else {
                int IG2 = IGNR[G + 1];
                L_CLR2 = (GFNR[G + 1][IG2] > GCMX[G + 1]);
            }

            double GCOWT;
            if (!L_CLR2) {
                // Upper layer member is cloudy
                if (!L_CLR1) {
                    GCOWT = FWTC[G][IG1];
                } else {
                    GCOWT = FWTC[G][0];
                }
            } else {
                // Upper layer member is clear
                if (!L_CLR1) {
                    GCOWT = FWTCC[G][IG1];
                } else {
                    GCOWT = FWTCC[G][0];
                }
            }
            WTCOL *= GCOWT;

            for (int L = GBOT[G]; L <= GTOP[G]; ++L) {
                if (NCLDF[L] >= GFNR[G][IG1]) {
                    ODCOL += CLT[L];
                }
            }
        }
        WCOL[I] = WTCOL;
        OCOL[I] = ODCOL;
    }
}

// =========================================================================
// ICA_QUD - Bin ICAs into NQD_ quadrature groups
// =========================================================================
/// Takes the full set of ICAs and groups them into NQD_ ranges of total OD.
/// Creates the Cumulative Probability Function and selects the mid-point ICA
/// for each group. The quadrature atmospheres have weights WTQCA.
///
/// Mode 6 = midpoints of each bin; Mode 7 = averaged properties within bins
/// (mode selection happens in the caller CLOUD_JX, not here).
///
/// Ported from ICA_QUD in src/Core/cldj_sub_mod.F90.
///
/// @param WCOL   Weight (fractional area) of each ICA [NICA used]
/// @param OCOL   Total cloud OD of each ICA [NICA used]
/// @param NICA   Number of ICAs
/// @param WTQCA  Output: quadrature weights [NQD_]
/// @param ISORT  Output: sorting index array [ICA_]
/// @param NQ1    Output: start index of each quadrature bin [NQD_]
/// @param NQ2    Output: end index of each quadrature bin [NQD_]
/// @param NDXQS  Output: mid-point index in each bin [NQD_]
inline void ICA_QUD(const double* WCOL, const double* OCOL,
                    int NICA,
                    double* WTQCA, int* ISORT,
                    int* NQ1, int* NQ2, int* NDXQS)
{
    // OD thresholds for quadrature bins (same as Fortran)
    constexpr double OD_QUAD[NQD_] = {0.5, 4.0, 30.0, 1.0e9};

    // Initialize outputs
    for (int i = 0; i < ICA_; ++i) ISORT[i] = 0;
    for (int i = 0; i < NQD_; ++i) {
        NQ1[i]   = 0;
        NQ2[i]   = 0;
        NDXQS[i] = 0;
        WTQCA[i] = 0.0;
    }

    // Sort ICAs by column OD (ascending)
    double OCOLS[ICA_];
    std::memset(OCOLS, 0, sizeof(OCOLS));

    if (NICA == 1) {
        ISORT[0] = 0;
        OCOLS[0] = OCOL[0];
    } else {
        HEAPSORT_A(NICA, OCOL, OCOLS, ISORT, ICA_);
    }

    // Build cumulative probability function
    double OCDFS[ICA_];
    OCDFS[0] = WCOL[ISORT[0]];
    for (int I = 1; I < NICA; ++I) {
        OCDFS[I] = OCDFS[I - 1] + WCOL[ISORT[I]];
    }

    // Find beginning/end of each quadrature range
    int I = 0;
    for (int N = 0; N < NQD_; ++N) {
        while (I < NICA && OCOLS[I] < OD_QUAD[N]) {
            ++I;
        }
        NQ2[N] = I - 1;
    }
    NQ1[0] = 0;
    for (int N = 1; N < NQD_; ++N) {
        NQ1[N] = NQ2[N - 1] + 1;
    }

    // Define QCA weights from cumulative prob, pick middle ICA
    for (int N = 0; N < NQD_; ++N) {
        int N1 = NQ1[N];
        int N2 = NQ2[N];
        if (N2 >= N1) {
            NDXQS[N] = (N1 + N2) / 2;
            if (N1 > 0) {
                WTQCA[N] = OCDFS[N2] - OCDFS[N1 - 1];
            } else {
                WTQCA[N] = OCDFS[N2];
            }
        }
    }
}

// =========================================================================
// ICA_III - Extract per-layer cloud OD for a specific ICA
// =========================================================================
inline void ICA_III(const double* CLT, int LTOP, int CBINU, int III,
                    const int* NCLDF,
                    const int GFNR[9][CBIN_ + 1],
                    const int* GNR, const int* GBOT, const int* GTOP,
                    int NRG, int NICA, double* TTCOL)
{
    (void)CBINU;  // used only for array bounds in Fortran
    std::memset(TTCOL, 0, LTOP * sizeof(double));

    int II = std::max(1, std::min(NICA, III));
    for (int G = 0; G < NRG; ++G) {
        int IG = (II - 1) % GNR[G];
        II = (II - 1) / GNR[G] + 1;
        for (int L = GBOT[G]; L <= GTOP[G]; ++L) {
            if (NCLDF[L] >= GFNR[G][IG]) {
                TTCOL[L] = CLT[L];
            }
        }
    }
}

} // namespace CloudOverlap
} // namespace CloudJ

// =========================================================================
// CLOUD_JX - Main cloud treatment driver
// =========================================================================
namespace CloudJ {

inline void CLOUD_JX(
    double U0, double SZA,
    const double RFL_flat[5 * (W_ + W_r)], double SOLF,
    bool LPRTJ,
    const double* PPP, const double* ZZZ,
    const double* TTT, const double* HHH,
    const double* DDD, const double* RRR,
    const double* OOO, const double* CCC,
    const double* LWP, const double* IWP,
    const double* REFFL, const double* REFFI,
    const double* CLDF, const int* CLDIW,
    double CLDCOR_in,
    const double* AERSP, const int* NDXAER,
    int L1U, int ANU, int NJXU,
    double* VALJXX, double* SKPERD, double* SWMSQ, double* OD18,
    int IRAN, int& NICA, int& JCOUNT, bool& LDARK,
    double* WTQCA,
    const CloudJState& state,
    int& rc,
    double* DirSfcFlux = nullptr,
    double* DiffSfcFlux = nullptr,
    double* DepFlux = nullptr,
    double* DiffTopFlux = nullptr)
{
    using namespace CloudOverlap;

    const int LWEPAR = state.LWEPAR;
    const int CLDFLAG = state.CLDFLAG;
    const int LNRG = state.LNRG;
    const int WW = W_ + W_r;

    rc = CLDJ_SUCCESS;
    bool LPRTJ0 = LPRTJ;
    JCOUNT = 0;
    NICA = 0;

    // Local working arrays
    std::vector<double> LWPX(L1U, 0.0), IWPX(L1U, 0.0);
    std::vector<double> REFFLX(L1U, 0.0), REFFIX(L1U, 0.0);

    // Zero primary outputs
    std::memset(VALJXX, 0, (L1U - 1) * NJXU * sizeof(double));
    std::memset(SKPERD, 0, (S_ + 2) * L1U * sizeof(double));
    std::memset(SWMSQ, 0, 6 * sizeof(double));
    std::memset(OD18, 0, L1U * sizeof(double));

    // Accumulation temporaries
    std::vector<double> VALJXXX((L1U - 1) * NJXU, 0.0);
    std::vector<double> SKPERDD((S_ + 2) * L1U, 0.0);
    double SWMSQQ[6] = {0};
    std::vector<double> OD18Q(L1U, 0.0);

    // Flux temporaries
    double FSBOT[W_ + W_r] = {0};
    double FJXBOT[W_ + W_r] = {0};
    std::vector<double> FLXD(L1U * WW, 0.0);
    std::vector<double> FJFLX(L1U * WW, 0.0);

    // Validate CLDFLAG
    if (CLDFLAG < 1 || CLDFLAG > 8 || CLDFLAG == 4) {
        CLOUDJ_ERROR(
            "Incorrect cloud index: must be 1-3 or 5-8",
            "CLOUD_JX in cloud_jx.hpp", rc);
        return;
    }

    // =====================================================================
    // CLDFLAG = 1, 2, 3: Simple cloud treatments, single PHOTO_JX call
    // =====================================================================
    if (CLDFLAG <= 3) {
        for (int L = 0; L < LWEPAR; ++L) {
            REFFLX[L] = REFFL[L];
            REFFIX[L] = REFFI[L];
        }

        if (CLDFLAG == 1) {
            // Clear sky - no clouds
            for (int L = 0; L < LWEPAR; ++L) {
                LWPX[L] = 0.0;
                IWPX[L] = 0.0;
                REFFLX[L] = 0.0;
                REFFIX[L] = 0.0;
            }
        } else if (CLDFLAG == 2) {
            // Average cloud cover: OD * CLDF
            for (int L = 0; L < LWEPAR; ++L) {
                double CLDFR = CLDF[L];
                LWPX[L] = LWP[L] * CLDFR;
                IWPX[L] = IWP[L] * CLDFR;
            }
        } else {
            // CLDFLAG == 3: OD * CLDF^1.5
            for (int L = 0; L < LWEPAR; ++L) {
                double CLDFR = CLDF[L] * std::sqrt(CLDF[L]);
                LWPX[L] = LWP[L] * CLDFR;
                IWPX[L] = IWP[L] * CLDFR;
            }
        }

        PhotoJX::PHOTO_JX(U0, SZA, RFL_flat, SOLF, LPRTJ0,
            PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
            LWPX.data(), IWPX.data(), REFFLX.data(), REFFIX.data(),
            AERSP, NDXAER,
            L1U, ANU, NJXU,
            VALJXX, SKPERD, SWMSQ, OD18,
            LDARK, FSBOT, FJXBOT, FLXD.data(), FJFLX.data(), state, rc);

        if (!LDARK) {
            JCOUNT = JCOUNT + 1;
        }

    } else {
        // =================================================================
        // CLDFLAG = 5, 6, 7, 8: ICA-based cloud overlap treatments
        // =================================================================

        // Prepare local cloud arrays
        std::vector<double> CLDX(LWEPAR, 0.0);
        std::vector<double> CLT(LWEPAR, 0.0);
        std::vector<double> CLTI(LWEPAR, 0.0);
        std::vector<double> CLTL(LWEPAR, 0.0);

        for (int L = 0; L < LWEPAR; ++L) {
            CLDX[L] = CLDF[L];
            LWPX[L] = LWP[L];
            IWPX[L] = IWP[L];
            REFFLX[L] = REFFL[L];
            REFFIX[L] = REFFI[L];
        }

        // Compute approximate cloud visible optical depths
        for (int L = 0; L < LWEPAR; ++L) {
            if (REFFIX[L] > 0.0) {
                CLTI[L] = IWPX[L] * 0.75 * 2.0 / (REFFIX[L] * 0.917);
                CLT[L] += CLTI[L];
            }
            if (REFFLX[L] > 0.0) {
                CLTL[L] = LWPX[L] * 0.75 * 2.1 / REFFLX[L];
                CLT[L] += CLTL[L];
            }
        }

        int LTOP = LWEPAR;

        // ICA overlap group setup
        std::vector<int> NCLDF(LTOP, 0);
        double CFBIN[CBIN_];
        int GFNR[9][CBIN_ + 1];
        int GCMX[9], GNR[9], GBOT[9], GTOP[9], GLVL[9];
        int NRG = 0;

        ICA_NR(CLDX.data(), CLT.data(), IWPX.data(), LWPX.data(),
               ZZZ, CLDIW,
               LTOP, CLDCOR_in, LNRG,
               NCLDF.data(), CFBIN, GFNR, GCMX, GNR, GBOT, GTOP, GLVL,
               NRG, NICA);

        // Generate all ICA weights and total ODs
        double WCOL[ICA_], OCOL[ICA_];
        ICA_ALL(CLDX.data(), CLT.data(), LTOP, CFBIN, CLDCOR_in,
                NCLDF.data(), GFNR, GCMX, GNR, GBOT, GTOP, GLVL,
                NRG, NICA, WCOL, OCOL);

        // Per-layer cloud OD for selected ICA
        std::vector<double> TTCOL(LTOP, 0.0);

        // ---- CLDFLAG = 5: Random ICA selection ----
        if (CLDFLAG == 5) {
            constexpr int NRANDO = 50;
            double WTRAN = 1.0 / static_cast<double>(NRANDO);

            // Build cumulative probability function
            double OCDFS[ICA_];
            OCDFS[0] = WCOL[0];
            for (int I = 1; I < NICA; ++I) {
                OCDFS[I] = OCDFS[I - 1] + WCOL[I];
            }

            for (int N = 0; N < NRANDO; ++N) {
                int IRANX = (IRAN + N) % NRAN_;
                double XRAN = static_cast<double>(state.RAN4[IRANX]);
                int I = 0;
                while (XRAN > OCDFS[I] && I < NICA - 1) {
                    ++I;
                }

                // Reconstruct per-layer cloud for ICA I+1 (1-based)
                ICA_III(CLT.data(), LTOP, CBIN_, I + 1,
                        NCLDF.data(), GFNR, GNR, GBOT, GTOP,
                        NRG, NICA, TTCOL.data());

                // Zero cloud paths not in selected ICA
                for (int L = 0; L < LTOP; ++L) {
                    LWPX[L] = LWP[L];
                    IWPX[L] = IWP[L];
                }
                for (int L = 0; L < LTOP; ++L) {
                    if (TTCOL[L] < 1.0e-8) {
                        LWPX[L] = 0.0;
                        IWPX[L] = 0.0;
                    }
                }

                PhotoJX::PHOTO_JX(U0, SZA, RFL_flat, SOLF, LPRTJ0,
                    PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
                    LWPX.data(), IWPX.data(), REFFLX.data(), REFFIX.data(),
                    AERSP, NDXAER,
                    L1U, ANU, NJXU,
                    VALJXXX.data(), SKPERDD.data(), SWMSQQ, OD18Q.data(),
                    LDARK, FSBOT, FJXBOT, FLXD.data(), FJFLX.data(),
                    state, rc);

                if (!LDARK) JCOUNT++;
                LPRTJ0 = false;

                // Accumulate weighted averages
                for (int J = 0; J < NJXU; ++J)
                    for (int L = 0; L < L1U - 1; ++L)
                        VALJXX[L + (L1U - 1) * J] +=
                            VALJXXX[L + (L1U - 1) * J] * WTRAN;
                for (int K = 0; K < S_ + 2; ++K)
                    for (int L = 0; L < L1U; ++L)
                        SKPERD[K + (S_ + 2) * L] +=
                            SKPERDD[K + (S_ + 2) * L] * WTRAN;
                for (int M = 0; M < 6; ++M)
                    SWMSQ[M] += SWMSQQ[M] * WTRAN;
                for (int L = 0; L < L1U; ++L)
                    OD18[L] += OD18Q[L] * WTRAN;
            }
        }

        // ---- CLDFLAG = 6: Quadrature midpoints ----
        if (CLDFLAG == 6) {
            int ISORT[ICA_], NQ1[NQD_], NQ2[NQD_], NDXQS[NQD_];
            ICA_QUD(WCOL, OCOL, NICA, WTQCA, ISORT, NQ1, NQ2, NDXQS);

            for (int N = 0; N < NQD_; ++N) {
                if (WTQCA[N] > 0.0) {
                    int I = ISORT[NDXQS[N]];

                    ICA_III(CLT.data(), LTOP, CBIN_, I + 1,
                            NCLDF.data(), GFNR, GNR, GBOT, GTOP,
                            NRG, NICA, TTCOL.data());

                    for (int L = 0; L < LTOP; ++L) {
                        LWPX[L] = LWP[L];
                        IWPX[L] = IWP[L];
                    }
                    for (int L = 0; L < LTOP; ++L) {
                        if (TTCOL[L] < 1.0e-8) {
                            LWPX[L] = 0.0;
                            IWPX[L] = 0.0;
                        }
                    }

                    PhotoJX::PHOTO_JX(U0, SZA, RFL_flat, SOLF, LPRTJ0,
                        PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
                        LWPX.data(), IWPX.data(), REFFLX.data(), REFFIX.data(),
                        AERSP, NDXAER,
                        L1U, ANU, NJXU,
                        VALJXXX.data(), SKPERDD.data(), SWMSQQ, OD18Q.data(),
                        LDARK, FSBOT, FJXBOT, FLXD.data(), FJFLX.data(),
                        state, rc);

                    if (!LDARK) JCOUNT++;
                    LPRTJ0 = false;

                    for (int J = 0; J < NJXU; ++J)
                        for (int L = 0; L < L1U - 1; ++L)
                            VALJXX[L + (L1U - 1) * J] +=
                                VALJXXX[L + (L1U - 1) * J] * WTQCA[N];
                    for (int K = 0; K < S_ + 2; ++K)
                        for (int L = 0; L < L1U; ++L)
                            SKPERD[K + (S_ + 2) * L] +=
                                SKPERDD[K + (S_ + 2) * L] * WTQCA[N];
                    for (int M = 0; M < 6; ++M)
                        SWMSQ[M] += SWMSQQ[M] * WTQCA[N];
                    for (int L = 0; L < L1U; ++L)
                        OD18[L] += OD18Q[L] * WTQCA[N];
                }
            }
        }

        // ---- CLDFLAG = 7: Quadrature averaged cloud ----
        if (CLDFLAG == 7) {
            int ISORT[ICA_], NQ1[NQD_], NQ2[NQD_], NDXQS[NQD_];
            ICA_QUD(WCOL, OCOL, NICA, WTQCA, ISORT, NQ1, NQ2, NDXQS);

            for (int N = 0; N < NQD_; ++N) {
                if (WTQCA[N] > 0.0 && NQ2[N] >= NQ1[N]) {
                    // Average cloud water paths over all ICAs in bin
                    std::fill(IWPX.begin(), IWPX.end(), 0.0);
                    std::fill(LWPX.begin(), LWPX.end(), 0.0);

                    for (int II = NQ1[N]; II <= NQ2[N]; ++II) {
                        int I = ISORT[II];
                        ICA_III(CLT.data(), LTOP, CBIN_, I + 1,
                                NCLDF.data(), GFNR, GNR, GBOT, GTOP,
                                NRG, NICA, TTCOL.data());

                        for (int L = 0; L < LTOP; ++L) {
                            if (TTCOL[L] > 1.0e-8) {
                                IWPX[L] += IWP[L] * WCOL[I];
                                LWPX[L] += LWP[L] * WCOL[I];
                            }
                        }
                    }

                    // Normalize by QCA weight
                    for (int L = 0; L < LTOP; ++L) {
                        IWPX[L] /= WTQCA[N];
                        LWPX[L] /= WTQCA[N];
                    }

                    PhotoJX::PHOTO_JX(U0, SZA, RFL_flat, SOLF, LPRTJ0,
                        PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
                        LWPX.data(), IWPX.data(), REFFLX.data(), REFFIX.data(),
                        AERSP, NDXAER,
                        L1U, ANU, NJXU,
                        VALJXXX.data(), SKPERDD.data(), SWMSQQ, OD18Q.data(),
                        LDARK, FSBOT, FJXBOT, FLXD.data(), FJFLX.data(),
                        state, rc);

                    if (!LDARK) JCOUNT++;
                    LPRTJ0 = false;

                    for (int J = 0; J < NJXU; ++J)
                        for (int L = 0; L < L1U - 1; ++L)
                            VALJXX[L + (L1U - 1) * J] +=
                                VALJXXX[L + (L1U - 1) * J] * WTQCA[N];
                    for (int K = 0; K < S_ + 2; ++K)
                        for (int L = 0; L < L1U; ++L)
                            SKPERD[K + (S_ + 2) * L] +=
                                SKPERDD[K + (S_ + 2) * L] * WTQCA[N];
                    for (int M = 0; M < 6; ++M)
                        SWMSQ[M] += SWMSQQ[M] * WTQCA[N];
                    for (int L = 0; L < L1U; ++L)
                        OD18[L] += OD18Q[L] * WTQCA[N];
                }
            }
        }

        // ---- CLDFLAG = 8: All ICAs ----
        if (CLDFLAG == 8) {
            for (int I = 0; I < NICA; ++I) {
                ICA_III(CLT.data(), LTOP, CBIN_, I + 1,
                        NCLDF.data(), GFNR, GNR, GBOT, GTOP,
                        NRG, NICA, TTCOL.data());

                for (int L = 0; L < LTOP; ++L) {
                    LWPX[L] = LWP[L];
                    IWPX[L] = IWP[L];
                }
                for (int L = 0; L < LTOP; ++L) {
                    if (TTCOL[L] < 1.0e-8) {
                        IWPX[L] = 0.0;
                        LWPX[L] = 0.0;
                    }
                }

                PhotoJX::PHOTO_JX(U0, SZA, RFL_flat, SOLF, LPRTJ0,
                    PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
                    LWPX.data(), IWPX.data(), REFFLX.data(), REFFIX.data(),
                    AERSP, NDXAER,
                    L1U, ANU, NJXU,
                    VALJXXX.data(), SKPERDD.data(), SWMSQQ, OD18Q.data(),
                    LDARK, FSBOT, FJXBOT, FLXD.data(), FJFLX.data(),
                    state, rc);

                if (!LDARK) JCOUNT++;
                LPRTJ0 = false;

                for (int J = 0; J < NJXU; ++J)
                    for (int L = 0; L < L1U - 1; ++L)
                        VALJXX[L + (L1U - 1) * J] +=
                            VALJXXX[L + (L1U - 1) * J] * WCOL[I];
                for (int K = 0; K < S_ + 2; ++K)
                    for (int L = 0; L < L1U; ++L)
                        SKPERD[K + (S_ + 2) * L] +=
                            SKPERDD[K + (S_ + 2) * L] * WCOL[I];
                for (int M = 0; M < 6; ++M)
                    SWMSQ[M] += SWMSQQ[M] * WCOL[I];
                for (int L = 0; L < L1U; ++L)
                    OD18[L] += OD18Q[L] * WCOL[I];
            }
        }
    }

    // Set optional flux outputs
    if (DirSfcFlux)
        std::memcpy(DirSfcFlux, FSBOT, WW * sizeof(double));
    if (DiffSfcFlux)
        std::memcpy(DiffSfcFlux, FJXBOT, WW * sizeof(double));
    if (DepFlux)
        std::memcpy(DepFlux, FLXD.data(), L1U * WW * sizeof(double));
    if (DiffTopFlux)
        std::memcpy(DiffTopFlux, FJFLX.data(), L1U * WW * sizeof(double));
}

} // namespace CloudJ

#endif // CLOUDJ_CLOUD_JX_HPP
