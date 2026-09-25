#ifndef CLOUDJ_STATE_HPP
#define CLOUDJ_STATE_HPP

#include <array>
#include <string>

namespace CloudJ {

// ===========================================================================
// Compile-time dimension parameters
// ===========================================================================

// Wavelength/spectral bin dimensions
constexpr int WX_ = 18;    // table dimension for cross-sections
constexpr int SX_ = 27;    // table dimension of broad-bands thru IR
constexpr int W_  = 18;    // no. of Fast-J Wavelength bins
constexpr int S_  = W_;    // v8.0: S_ == W_ == 18 (no Solar-J extension)
constexpr int NW1 = 1;
constexpr int NW2 = W_;
constexpr int NS1 = 1;
constexpr int NS2 = S_;
constexpr int W_r = S_ - W_;  // = 0 for v8.0

// Cross-section and aerosol dimensions (model-dependent)
#ifdef MODEL_GEOSCHEM
constexpr int X_  = 123;   // max no. of X-section data sets
constexpr int A_  = 56;    // max no. of Aerosol Mie sets
constexpr int AN_ = 37;    // max no. of FJX aerosols in layer
#elif defined(MODEL_STANDALONE)
constexpr int X_  = 72;
constexpr int A_  = 40;
constexpr int AN_ = 25;
#else
constexpr int X_  = 72;
constexpr int A_  = 40;
constexpr int AN_ = 25;
#endif

// Cloud data dimensions
constexpr int C_  = 3;     // no. of cloud data sets (liquid, irreg-ice, hex-ice)
constexpr int CR_ = 6;     // no. of effective radii per cloud data set

// Mie scattering array levels
constexpr int N_  = 601;   // no. of levels in Mie scattering arrays

// Gauss quadrature
constexpr int M_  = 4;     // no. of Gauss points (8-stream)
constexpr int M2_ = 2 * M_;  // = 8, replaces MFIT

// J-value dimension
constexpr int JVN_ = 200;  // max no. of J-values

// Stratospheric sulfate aerosol / GeoMIP dimensions
constexpr int SSA_DIM = 18;   // SSA_ in Fortran
constexpr int GGA_    = 15;

// Cloud overlap parameters
constexpr int CBIN_ = 10;    // no. of quantized cloud fraction bins
constexpr int ICA_  = 20000; // max no. of independent column atmospheres
constexpr int NQD_  = 4;     // no. of cloud quadrature bins

// Random number dimension
constexpr int NRAN_ = 10007;

// Reference profile dimensions
constexpr int LREF  = 51;    // layer dim in reference profiles
constexpr int JREF  = 18;    // latitude dim in reference profiles
constexpr int LGREF = 19;    // layer dim in GeoMIP reference profiles

// Spectral sub-bin data (NGC)
constexpr int NSBIN = 27;

// ===========================================================================
// Physical constants (identical to Fortran reference values)
// ===========================================================================

constexpr double RAD      = 6375.0e5;     // Radius of Earth (cm)
constexpr double ZZHT     = 5.0e5;        // Scale height (cm) above top of CTM
constexpr double MASFAC   = 100.0 * 6.022e+23 / (28.97 * 9.8 * 10.0);  // pressure->column density
constexpr double HeatFac_ = 86400.0 * 9.80616 / 1.00464e5;  // W/m2 -> K/day
constexpr double CPI      = 3.141592653589793;
constexpr double C2PI     = 2.0 * CPI;
constexpr double CPI180   = CPI / 180.0;
constexpr double G0       = 9.80665;      // standard gravity (m/s2)
constexpr double G100     = 100.0 / G0;

// ===========================================================================
// 4-point Gauss quadrature arrays (8-stream)
// ===========================================================================

constexpr std::array<double, M_> EMU = {
    0.06943184420297, 0.33000947820757,
    0.66999052179243, 0.93056815579703
};

constexpr std::array<double, M_> WT = {
    0.17392742256873, 0.32607257743127,
    0.32607257743127, 0.17392742256873
};

// ===========================================================================
// NGC: sub-bin counts per wavelength bin (all 1 for Cloud-J v8.0)
// ===========================================================================

constexpr std::array<int, 27> NGC = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1
};

// ===========================================================================
// Cloud flag descriptions
// ===========================================================================

constexpr int NUM_CLDFLAG_OPTIONS = 8;

// ===========================================================================
// CloudJState struct — encapsulates all shared data from cldj_cmn_mod.F90
// ===========================================================================

struct CloudJState {
    // --- Grid dimensions ---
#ifdef MODEL_STANDALONE
    static constexpr int L_      = 57;   // no. of CTM layers
    static constexpr int L1_     = 58;   // L_+1
    static constexpr int L2_     = 59;   // L_+2
    static constexpr int LWEPAR  = 34;   // no. layers with clouds
#else
    int L_      = 0;    // no. of CTM layers (set at runtime)
    int L1_     = 0;    // L_+1
    int L2_     = 0;    // L_+2
    int LWEPAR  = 0;    // no. layers with clouds
#endif

    // --- Runtime configuration ---
    int    ATM0     = 0;       // spherical correction option (0=flat,1=sphr,2=refr,3=geom)
    double ATAU     = 1.05;    // factor increase in cloud OD layer-to-layer
    double ATAU0    = 0.005;   // minimum cloud OD in uppermost inserted layer
    bool   USEH2OUV = false;   // whether to use H2O UV absorption
    int    CLDFLAG  = 7;       // cloud treatment flag (1-8)
    double CLDCOR   = 0.33;    // cloud decorrelation between max-overlap blocks
    int    LNRG     = 6;       // number of max-overlap blocks
    int    NWBIN    = 18;      // number of wavelength bins to use
    int    NJX      = 0;       // no. of fast-JX J-values (set from cross-section data)

    // =========================================================================
    // Spectral data (from FJX_spec.dat via RD_XXX)
    // =========================================================================

    // WL: Centres of wavelength bins - effective wavelength (nm)
    std::array<double, SX_> WL{};

    // WBIN: Boundaries of wavelength bins (microns)
    std::array<double, SX_ + 1> WBIN{};

    // FL: Solar flux incident on top of atmosphere (cm-2.s-1)
    std::array<double, SX_> FL{};

    // FW: Solar flux in W/m2
    std::array<double, SX_> FW{};

    // FPAR: PAR quantum action spectrum
    std::array<double, SX_> FPAR{};

    // QRAYL: Rayleigh parameters - effective cross-section (cm2)
    std::array<double, SX_> QRAYL{};

    // SJSUB: breakdown of super-bins into sub-bins
    double SJSUB[SX_][16]{};

    // QH2O: H2O UV-blue cross-sections (290-350 nm)
    std::array<double, WX_> QH2O{};

    // KDOKR, LDOKR: RRTMG index arrays
    int KDOKR[100]{};
    int LDOKR[100]{};

    // NSJSUB: sub-bin counts
    int NSJSUB[SX_]{};

    // QO2: O2 cross-sections [3][WX_] — first index is temperature node, second is wavelength bin
    double QO2[3][WX_]{};

    // QO3: O3 cross-sections [3][WX_] — first index is temperature node, second is wavelength bin
    double QO3[3][WX_]{};

    // Q1D: O3 => O(1D) quantum yield [3][WX_] — first index is temperature node, second is wavelength bin
    double Q1D[3][WX_]{};

    // QQQ: Supplied cross sections in each wavelength bin (cm2) [WX_][3][X_]
    double QQQ[WX_][3][X_]{};

    // TQQ: Temperature for supplied cross sections [3][X_]
    double TQQ[3][X_]{};

    // LQQ: interpolation flag (1,2,3 = T or P) [X_]
    int LQQ[X_]{};

    // SQQ: Flag (pressure or temperature tables) [X_]
    char SQQ[X_]{};

    // TITLEJX: Short title for supplied cross sections
    std::array<std::string, X_> TITLEJX{};

    // TITLEJL: Long title for supplied cross sections
    std::array<std::string, X_> TITLEJL{};

    // =========================================================================
    // Aerosol Mie data (from FJX_scat-aer.dat via RD_MIE)
    // =========================================================================

    // NAA: Number of aerosol scattering categories
    int NAA = 0;

    // TITLAA: Aerosol Mie Titles
    std::array<std::string, A_> TITLAA{};

    // QAA: Aerosol Q-ext [5][A_]
    double QAA[5][A_]{};

    // WAA: Wavelengths for phase functions [5][A_]
    double WAA[5][A_]{};

    // PAA: Phase function expansion (8 terms) [8][5][A_]
    double PAA[8][5][A_]{};

    // RAA: Effective radius per aerosol type [A_]
    double RAA[A_]{};

    // SAA: Single scattering albedo [5][A_]
    double SAA[5][A_]{};

    // DAA: Density (g/cm3) [A_]
    double DAA[A_]{};

    // =========================================================================
    // Cloud scattering data (from FJX_scat-cld.dat via RD_CLD)
    // =========================================================================

    // NCC, MCC: Number of cloud categories and effective radii
    int NCC = 0;
    int MCC = 0;

    // TITLCC: Cloud type titles
    std::array<std::string, C_> TITLCC{};

    // RCC: Effective radius per cloud type [CR_][C_]
    double RCC[CR_][C_]{};

    // GCC: Effective geometric cross section [CR_][C_]
    double GCC[CR_][C_]{};

    // DCC: Density (g/cm3) [C_]
    double DCC[C_]{};

    // QCC: Cloud Q-ext [SX_][CR_][C_]
    double QCC[SX_][CR_][C_]{};

    // WCC: Wavelengths for cloud phase functions [SX_][C_]
    double WCC[SX_][C_]{};

    // SCC: Cloud single scattering albedo [SX_][CR_][C_]
    double SCC[SX_][CR_][C_]{};

    // PCC: Cloud phase function expansion [8][SX_][CR_][C_]
    double PCC[8][SX_][CR_][C_]{};

    // =========================================================================
    // Stratospheric Sulfate Aerosol data (from FJX_scat-ssa.dat via RD_SSA)
    // =========================================================================

    // NSS: Number of SSA categories
    int NSS = 0;

    // TITLSS: SSA type titles
    std::array<std::string, SSA_DIM> TITLSS{};

    // RSS: Effective radius [SSA_DIM]
    double RSS[SSA_DIM]{};

    // GSS: Effective geometric cross section [SSA_DIM]
    double GSS[SSA_DIM]{};

    // DSS: Density (g/cm3) [SSA_DIM]
    double DSS[SSA_DIM]{};

    // TSS: Temperature (K) [SSA_DIM]
    double TSS[SSA_DIM]{};

    // WSS: Weight percent sulfuric acid [SSA_DIM]
    double WSS[SSA_DIM]{};

    // QSS: Q-ext [SX_][SSA_DIM]
    double QSS[SX_][SSA_DIM]{};

    // SSS: Single scattering albedo [SX_][SSA_DIM]
    double SSS[SX_][SSA_DIM]{};

    // PSS: Phase function expansion [8][SX_][SSA_DIM]
    double PSS[8][SX_][SSA_DIM]{};

    // =========================================================================
    // GeoMIP scattering data (from FJX_scat-geo.dat via RD_GEO)
    // =========================================================================

    // NGG: Number of GeoMIP categories
    int NGG = 0;

    // RGG: Effective radius [GGA_]
    double RGG[GGA_]{};

    // DGG: Density (g/cm3) [GGA_]
    double DGG[GGA_]{};

    // QGG: Q-ext [SX_][GGA_]
    double QGG[SX_][GGA_]{};

    // SGG: Single scattering albedo [SX_][GGA_]
    double SGG[SX_][GGA_]{};

    // PGG: Phase function expansion [8][SX_][GGA_]
    double PGG[8][SX_][GGA_]{};

    // =========================================================================
    // UMich aerosol data (from FJX_scat-UMa.dat via RD_UM)
    // =========================================================================

    // WMM: U Michigan aerosol wavelengths [6]
    double WMM[6]{};

    // UMAER: U Michigan aerosol data sets [3][6][21][33]
    double UMAER[3][6][21][33]{};

    // =========================================================================
    // Climatology reference profiles (from atmos_std.dat, atmos_h2och4.dat)
    // NOTE: only used in Cloud-J standalone
    // =========================================================================

    // T_REF: Temperature reference profiles [LREF][JREF][12]
    double T_REF[LREF][JREF][12]{};

    // O_REF: Ozone reference profiles [LREF][JREF][12]
    double O_REF[LREF][JREF][12]{};

    // H2O_REF: Water vapor reference profiles [LREF][JREF][12]
    double H2O_REF[LREF][JREF][12]{};

    // CH4_REF: Methane reference profiles [LREF][JREF][12]
    double CH4_REF[LREF][JREF][12]{};

    // =========================================================================
    // GeoMIP reference profiles (from atmos_geomip.dat)
    // =========================================================================

    // R_GREF: Reff (microns) [64][LGREF][12]
    double R_GREF[64][LGREF][12]{};

    // X_GREF: micro-g-H2SO4/kg-air [64][LGREF][12]
    double X_GREF[64][LGREF][12]{};

    // A_GREF: 4*pi*R^2 = microns^2/cm^3 [64][LGREF][12]
    double A_GREF[64][LGREF][12]{};

    // Y_GREF: latitude array [64]
    double Y_GREF[64]{};

    // P_GREF: pressure levels [LGREF]
    double P_GREF[LGREF]{};

    // =========================================================================
    // J-value mapping data (from FJX_j2j.dat via RD_JS_JX)
    // =========================================================================

    // JFACTA: multiplication factor for calculated J [JVN_]
    double JFACTA[JVN_]{};

    // JIND: index arrays mapping Jvalue(j) onto rates [JVN_]
    int JIND[JVN_]{};

    // NRATJ: number of photolysis reactions in CTM chemistry
    int NRATJ = 0;

    // JVMAP: label of J-value used to match with FJX J's
    std::array<std::string, JVN_> JVMAP{};

    // JLABEL: label of J-value used in the chem model
    std::array<std::string, JVN_> JLABEL{};

    // BRANCH: Branches for photolysis species [JVN_]
    int BRANCH[JVN_]{};

    // RNAMES: Names of photolysis species
    std::array<std::string, JVN_> RNAMES{};

    // =========================================================================
    // Random number array for cloud-JX
    // =========================================================================

    // RAN4: Random number set (single precision, matching Fortran real*4)
    std::array<float, NRAN_> RAN4{};
};

} // namespace CloudJ

#endif // CLOUDJ_STATE_HPP
