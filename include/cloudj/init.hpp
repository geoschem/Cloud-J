#ifndef CLOUDJ_INIT_HPP
#define CLOUDJ_INIT_HPP

#include <cloudj/state.hpp>
#include <cloudj/error.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>

// Table data headers
#include <cloudj/tables/FJX_spec.hpp>
#include <cloudj/tables/FJX_scat_cld.hpp>
#include <cloudj/tables/FJX_scat_ssa.hpp>
#include <cloudj/tables/FJX_scat_aer.hpp>
#include <cloudj/tables/FJX_scat_UMa.hpp>
#include <cloudj/tables/FJX_scat_geo.hpp>
#include <cloudj/tables/FJX_j2j.hpp>
#include <cloudj/tables/atmos_std.hpp>
#include <cloudj/tables/atmos_h2och4.hpp>
#include <cloudj/tables/atmos_geomip.hpp>

namespace CloudJ::Init {

// =========================================================================
// load_table_or_default: on-disk table override with embedded-table fallback
// =========================================================================
//
// Returns the contents of `datadir/filename` if datadir is non-empty and the
// file exists, opens, and is non-empty; otherwise returns `embedded_default`
// unchanged. Prints one line to stderr only when a disk file is actually
// used, so silent embedded-table use (the default) produces no output.
inline std::string load_table_or_default(const std::string& datadir,
                                          const std::string& filename,
                                          const char* embedded_default) {
    if (!datadir.empty()) {
        std::string path = datadir;
        if (!path.empty() && path.back() != '/') path += '/';
        path += filename;
        std::ifstream f(path, std::ios::in | std::ios::binary);
        if (f.is_open()) {
            std::string contents((std::istreambuf_iterator<char>(f)),
                                 std::istreambuf_iterator<char>());
            if (!contents.empty()) {
                std::cerr << "Cloud-J: loaded " << filename << " from "
                          << datadir << "\n";
                return contents;
            }
        }
    }
    return std::string(embedded_default);
}

// =========================================================================
// RANSET: generates pseudo-random numbers using Knuth's subtractive method
// Ported from cldj_init_mod.F90 RANSET subroutine
// =========================================================================

inline void RANSET(std::array<float, NRAN_>& RAN4, int& ISTART, int& rc) {
    constexpr int MBIG = 1000000000;
    constexpr int MSEED = 161803398;
    constexpr int MZ = 0;
    constexpr float FAC = 1.0e-9f;

    rc = CLDJ_SUCCESS;

    int MA[55];
    int MJ = MSEED - std::abs(ISTART);
    MJ = MJ % MBIG;
    MA[54] = MJ;  // MA(55) in Fortran -> MA[54] in C++
    int MK = 1;

    for (int I = 1; I <= 54; ++I) {
        int II = (21 * I) % 55;
        MA[II - 1] = MK;  // Fortran 1-indexed
        MK = MJ - MK;
        if (MK < MZ) MK += MBIG;
        MJ = MA[II - 1];
    }
    for (int K = 1; K <= 4; ++K) {
        for (int I = 1; I <= 55; ++I) {
            MA[I - 1] = MA[I - 1] - MA[((I + 30) % 55)];
            if (MA[I - 1] < MZ) MA[I - 1] += MBIG;
        }
    }
    int INEXT = 0;
    int INEXTP = 31;
    ISTART = 1;

    // Generate NRAN_ pseudo-random numbers
    for (int J = 0; J < NRAN_; ++J) {
        INEXT = (INEXT % 55) + 1;
        INEXTP = (INEXTP % 55) + 1;
        MJ = MA[INEXT - 1] - MA[INEXTP - 1];
        if (MJ < MZ) MJ += MBIG;
        MA[INEXT - 1] = MJ;
        RAN4[J] = static_cast<float>(MJ) * FAC;
    }
}

// =========================================================================
// RD_XXX: Read spectral data from FJX_spec table
// Ported from cldj_init_mod.F90
// =========================================================================

inline void RD_XXX(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_XXX (C++ port)";
    rc = CLDJ_SUCCESS;

    // Initialize TQQ
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < X_; ++j)
            state.TQQ[i][j] = 0.0;

    if (W_ != 18) {
        CLOUDJ_ERROR("no. wavelengths wrong: W_ != 18", thisloc, rc);
        return;
    }

    std::istringstream iss(load_table_or_default(datadir, "FJX_spec.dat", Tables::FJX_spec));
    std::string line;

    // Line 1: title
    std::getline(iss, line);  // TIT_SPEC

    // Line 2: blank/comment
    std::getline(iss, line);

    // Line 3: NWWW, NSSS
    std::getline(iss, line);
    int NWWW = 0, NSSS = 0;
    std::istringstream dim_ss(line);
    dim_ss >> NWWW;
    // skip text until next number
    std::string tmp;
    while (dim_ss >> tmp) {
        try { NSSS = std::stoi(tmp); break; }
        catch (...) { continue; }
    }

    if (NWWW != WX_ || NSSS != SX_) {
        CLOUDJ_ERROR("WX_ or SX_ incompatible data sets", thisloc, rc);
        return;
    }

    // Helper lambda to read data values from lines with 5-char prefix
    // Reads NSSS values spread across multiple lines (6 per line after 5-char prefix)
    auto read_spectral_array = [&](double* arr, int count) {
        int read = 0;
        while (read < count) {
            std::getline(iss, line);
            // Skip the 5-character prefix (e.g., "01:06" or "    a")
            std::istringstream ls(line.substr(5));
            double val;
            while (ls >> val && read < count) {
                arr[read++] = val;
            }
        }
    };

    // Read w-eff (WL) - header line then data
    std::getline(iss, line);  // header: "w-eff !wavelength..."
    read_spectral_array(state.WL.data(), NSSS);

    // Read w-bins (WBIN) - header line then data
    std::getline(iss, line);  // header: "w-bins!wavel..."
    read_spectral_array(state.WBIN.data(), NSSS);

    // Read SPhot (FL) - header line then data
    std::getline(iss, line);  // header: "SPhot !solar..."
    read_spectral_array(state.FL.data(), NSSS);

    // Read SWatt (FW) - header line then data
    std::getline(iss, line);  // header: "SWatt |solarheat..."
    read_spectral_array(state.FW.data(), NSSS);

    // Read Y-PAR (FPAR) - header line then data
    std::getline(iss, line);  // header: "Y-PAR !photosyn..."
    read_spectral_array(state.FPAR.data(), NSSS);

    // Read X-Rayl (QRAYL) - header line then data
    std::getline(iss, line);  // header: "X-Rayl!Rayleigh..."
    read_spectral_array(state.QRAYL.data(), NSSS);

    // Read SJ sub-bins (SJSUB) - header line then data
    // In Cloud-J v8.0, S_==W_==18, so these are read but then reset to 1/0
    std::getline(iss, line);  // header: "SJbins|SJ bin fractions..."
    for (int I = NWWW - 1; I < NSSS; ++I) {  // Fortran: do I=NWWW, NSSS
        // Each sub-bin has 3 lines of 5 values each (15 total)
        int sub_read = 0;
        for (int row = 0; row < 3; ++row) {
            std::getline(iss, line);
            std::istringstream ls(line.substr(5));
            double val;
            while (ls >> val && sub_read < 15) {
                state.SJSUB[I][sub_read++] = val;
            }
        }
    }
    // Reset SJSUB for Cloud-J (only first sub-bin used)
    for (int I = 0; I < NSSS; ++I) {
        state.SJSUB[I][0] = 1.0;
        for (int IW = 1; IW < 15; ++IW) {
            state.SJSUB[I][IW] = 0.0;
        }
    }

    // Read H2O cross-sections (QH2O)
    std::getline(iss, line);  // header: "H2O   !H2O UV absorpt..."
    {
        int read = 0;
        while (read < NWWW) {
            std::getline(iss, line);
            std::istringstream ls(line.substr(5));
            double val;
            while (ls >> val && read < NWWW) {
                state.QH2O[read++] = val;
            }
        }
    }
    // If H2O UV absorption turned off, zero out QH2O
    if (!state.USEH2OUV) {
        state.QH2O.fill(0.0);
    }

    // Helper to read cross-section block: SQQ char, TQQ temp, then 18 values over 3 lines
    auto read_xsect_block = [&](char& sqq, double& tqq, double* xsect) {
        std::getline(iss, line);
        // Format: "x300a val1 val2 val3 val4 val5 val6"
        // First char is SQQ flag, next 3 chars are temperature, then 'a'/'b'/'c'
        sqq = line[0];
        tqq = std::stod(line.substr(1, 3));
        std::istringstream ls(line.substr(5));
        int read = 0;
        double val;
        while (ls >> val && read < 6) {
            xsect[read++] = val;
        }
        // Second line (6 values)
        std::getline(iss, line);
        std::istringstream ls2(line.substr(5));
        while (ls2 >> val && read < 12) {
            xsect[read++] = val;
        }
        // Third line (6 values)
        std::getline(iss, line);
        std::istringstream ls3(line.substr(5));
        while (ls3 >> val && read < 18) {
            xsect[read++] = val;
        }
    };

    // Helper to read continuation block (no SQQ, just temp + data)
    auto read_xsect_cont = [&](double& tqq, double* xsect) {
        std::getline(iss, line);
        // Format: " 260a val1 val2..."
        tqq = std::stod(line.substr(1, 3));
        std::istringstream ls(line.substr(5));
        int read = 0;
        double val;
        while (ls >> val && read < 6) {
            xsect[read++] = val;
        }
        std::getline(iss, line);
        std::istringstream ls2(line.substr(5));
        while (ls2 >> val && read < 12) {
            xsect[read++] = val;
        }
        std::getline(iss, line);
        std::istringstream ls3(line.substr(5));
        while (ls3 >> val && read < 18) {
            xsect[read++] = val;
        }
    };

    // Read O2 cross-sections (3 temperatures): header + 3 blocks
    {
        std::getline(iss, line);  // header "O2    !O2=O+O..."
        char sqq_tmp;
        read_xsect_block(sqq_tmp, state.TQQ[0][0], state.QO2[0]);
        state.SQQ[0] = sqq_tmp;
        // Read continuation for T2
        std::getline(iss, line);  // "O2" repeat title
        read_xsect_cont(state.TQQ[1][0], state.QO2[1]);
        // Read continuation for T3
        std::getline(iss, line);  // "O2" repeat title
        read_xsect_cont(state.TQQ[2][0], state.QO2[2]);
        state.TITLEJX[0] = "O2    ";
        state.TITLEJL[0] = "O2=O+O          ";
        state.LQQ[0] = 3;
    }

    // Read O3 cross-sections (3 temperatures)
    {
        std::getline(iss, line);  // header "O3    !O3-total..."
        char sqq_tmp;
        read_xsect_block(sqq_tmp, state.TQQ[0][1], state.QO3[0]);
        state.SQQ[1] = sqq_tmp;
        std::getline(iss, line);  // "O3" repeat
        read_xsect_cont(state.TQQ[1][1], state.QO3[1]);
        std::getline(iss, line);  // "O3" repeat
        read_xsect_cont(state.TQQ[2][1], state.QO3[2]);
        state.TITLEJX[1] = "O3    ";
        state.TITLEJL[1] = "O3-total        ";
        state.LQQ[1] = 3;
    }

    // Read O3(1D) quantum yields (3 temperatures)
    {
        std::getline(iss, line);  // header "O3(1D)!Qyld..."
        char sqq_tmp;
        read_xsect_block(sqq_tmp, state.TQQ[0][2], state.Q1D[0]);
        state.SQQ[2] = sqq_tmp;
        std::getline(iss, line);  // "O3(1D)" repeat
        read_xsect_cont(state.TQQ[1][2], state.Q1D[1]);
        std::getline(iss, line);  // "O3(1D)" repeat
        read_xsect_cont(state.TQQ[2][2], state.Q1D[2]);
        state.TITLEJX[2] = "O3(1D)";
        state.TITLEJL[2] = "Qyld O3=O(1D)+O2";
        state.LQQ[2] = 3;
    }

    // Read remaining species: variable number of T/P sets per J-value
    int JJ = 3;  // Already read O2, O3, O3(1D) as indices 0,1,2

    // Read the first header for the next species
    std::getline(iss, line);  // Next species header
    std::string TIT_J1S = line.substr(0, 6);
    std::string TIT_J1L = (line.size() > 7) ? line.substr(7, 16) : "";

    while (TIT_J1S != "endofJ") {
        if (JJ >= X_) {
            CLOUDJ_ERROR("X_ not large enough", thisloc, rc);
            return;
        }

        state.TITLEJX[JJ] = TIT_J1S;
        state.TITLEJL[JJ] = TIT_J1L;

        // Read first T/P block
        char sqq_tmp;
        double tqq_val;
        double xsect[WX_] = {};
        read_xsect_block(sqq_tmp, tqq_val, xsect);
        state.SQQ[JJ] = sqq_tmp;
        state.TQQ[0][JJ] = tqq_val;
        for (int IW = 0; IW < NWWW; ++IW)
            state.QQQ[IW][0][JJ] = xsect[IW];
        state.LQQ[JJ] = 1;

        // Try to read next header (might be same species for 2nd T/P)
        std::getline(iss, line);
        if (line.empty() || iss.eof()) break;
        std::string next_title = line.substr(0, 6);
        std::string next_long = (line.size() > 7) ? line.substr(7, 16) : "";

        if (next_title == "endofJ") break;

        if (next_title == state.TITLEJX[JJ]) {
            // Read 2nd T/P
            read_xsect_cont(tqq_val, xsect);
            state.TQQ[1][JJ] = tqq_val;
            for (int IW = 0; IW < NWWW; ++IW)
                state.QQQ[IW][1][JJ] = xsect[IW];
            state.LQQ[JJ] = 2;

            // Try 3rd T/P
            std::getline(iss, line);
            if (line.empty() || iss.eof()) { JJ++; break; }
            next_title = line.substr(0, 6);
            next_long = (line.size() > 7) ? line.substr(7, 16) : "";

            if (next_title == "endofJ") { JJ++; break; }

            if (next_title == state.TITLEJX[JJ]) {
                read_xsect_cont(tqq_val, xsect);
                state.TQQ[2][JJ] = tqq_val;
                for (int IW = 0; IW < NWWW; ++IW)
                    state.QQQ[IW][2][JJ] = xsect[IW];
                state.LQQ[JJ] = 3;

                // Read next header
                std::getline(iss, line);
                if (line.empty() || iss.eof()) { JJ++; break; }
                TIT_J1S = line.substr(0, 6);
                TIT_J1L = (line.size() > 7) ? line.substr(7, 16) : "";
            } else {
                TIT_J1S = next_title;
                TIT_J1L = next_long;
            }
        } else {
            TIT_J1S = next_title;
            TIT_J1L = next_long;
        }
        JJ++;
    }

    state.NJX = JJ;

    // TROP-ONLY: drop strat X-sections if NWBIN < 18
    if (state.NWBIN == 12 || state.NWBIN == 8) {
        int jj_new = 3;
        for (int J = 3; J < state.NJX; ++J) {
            if (state.SQQ[J] != 'x') {
                if (jj_new < J) {
                    state.TITLEJX[jj_new] = state.TITLEJX[J];
                    state.TITLEJL[jj_new] = state.TITLEJL[J];
                    state.LQQ[jj_new] = state.LQQ[J];
                    state.SQQ[jj_new] = state.SQQ[J];
                    for (int LQ = 0; LQ < state.LQQ[J]; ++LQ) {
                        state.TQQ[LQ][jj_new] = state.TQQ[LQ][J];
                        for (int IW = 0; IW < NWWW; ++IW)
                            state.QQQ[IW][LQ][jj_new] = state.QQQ[IW][LQ][J];
                    }
                }
                jj_new++;
            }
        }
        state.NJX = jj_new;
    }

    // Need to check that TQQ (= T(K) or p(hPa)) is monotonically increasing
    // (Fortran RD_XXX in cldj_init_mod.F90, after the TROP-only collapse).
    for (int J = 0; J < state.NJX; ++J) {
        if ((state.LQQ[J] == 3) && (state.TQQ[1][J] >= state.TQQ[2][J])) {
            CLOUDJ_ERROR("TQQ out of order", thisloc, rc);
            return;
        }
        if ((state.LQQ[J] == 2) && (state.TQQ[0][J] >= state.TQQ[1][J])) {
            CLOUDJ_ERROR("TQQ out of order", thisloc, rc);
            return;
        }
    }

    // Zero FL for reduced wavelength bins
    if (state.NWBIN == 12) {
        for (int IW = 0; IW < 4; ++IW) state.FL[IW] = 0.0;
        state.FL[8] = 0.0;
        state.FL[9] = 0.0;
    }
    if (state.NWBIN == 8) {
        for (int IW = 0; IW < 4; ++IW) state.FL[IW] = 0.0;
        state.FL[4] *= 2.0;
        for (int IW = 5; IW < 11; ++IW) state.FL[IW] = 0.0;
    }
}

// =========================================================================
// RD_CLD: Read cloud scattering data from FJX_scat_cld table
// =========================================================================

inline void RD_CLD(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_CLD (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "FJX_scat-cld.dat", Tables::FJX_scat_cld));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2: NCC
    std::getline(iss, line);
    state.NCC = std::stoi(line);
    // Line 3: MCC
    std::getline(iss, line);
    state.MCC = std::stoi(line);
    // Skip 5 comment lines
    for (int i = 0; i < 5; ++i) std::getline(iss, line);

    for (int K = 0; K < state.NCC; ++K) {
        // Read cloud type header: title + density
        std::getline(iss, line);
        state.TITLCC[K] = line.substr(0, 12);
        state.DCC[K] = std::stod(line.substr(12, 8));

        // Read data for each wavelength bin from 12 to SX_ (Fortran: J=12,SX_)
        // The reference reads each row with a fixed-column format
        //   (i2,1x,f5.2,f5.1,f7.1,f5.3,e8.1,f6.3,f8.5,7f6.3)
        // so fields must be sliced by column, not split on whitespace. Some
        // rows (e.g. the largest ice effective radius) jam two numeric fields
        // together with no separator ("95.519927." = Reff 95.5 + cross-section
        // 19927.); a token-based read mis-parses these and silently zeroes the
        // trailing phase-function values.
        for (int J = 11; J < SX_; ++J) {
            for (int I = 0; I < state.MCC; ++I) {
                std::getline(iss, line);
                // Pad so short/blank-padded rows still yield every field.
                line.resize(90, ' ');
                auto fld = [&line](int pos, int width) -> double {
                    return std::stod(line.substr(pos, width));
                };
                int JCC = static_cast<int>(fld(0, 2));
                if (JCC != J + 1) {
                    // Matches the Fortran err/stop path on a malformed row.
                    CLOUDJ_ERROR("Error in read", thisloc, rc);
                    return;
                }
                state.WCC[J][K] = fld(3, 5);
                state.RCC[I][K] = fld(8, 5);
                state.GCC[I][K] = fld(13, 7);
                // Reff index (XNDR) and imaginary index (XNDI) are not used
                // downstream, but occupy columns 20-24 and 25-32.
                state.QCC[J][I][K] = fld(33, 6);
                state.SCC[J][I][K] = fld(39, 8);
                // Columns 47 onward: 7 values -> PCC[1..7] (Fortran L=2..8).
                for (int L = 1; L < 8; ++L) {
                    state.PCC[L][J][I][K] = fld(47 + 6 * (L - 1), 6);
                }
                state.PCC[0][J][I][K] = 1.0;
            }
            // Skip blank line between wavelength groups
            std::getline(iss, line);
        }
    }

    // Replicate cloud data for wavelengths < 295 nm (J=0..10) from J=11
    for (int K = 0; K < state.NCC; ++K) {
        for (int J = 0; J < 11; ++J) {
            state.WCC[J][K] = state.WCC[11][K];
            for (int I = 0; I < state.MCC; ++I) {
                state.QCC[J][I][K] = state.QCC[11][I][K];
                state.SCC[J][I][K] = state.SCC[11][I][K];
                for (int L = 0; L < 8; ++L) {
                    state.PCC[L][J][I][K] = state.PCC[L][11][I][K];
                }
            }
        }
    }
}

// =========================================================================
// RD_SSA: Read stratospheric sulfate aerosol scattering data
// =========================================================================

inline void RD_SSA(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_SSA (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "FJX_scat-ssa.dat", Tables::FJX_scat_ssa));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2: blank/comment
    std::getline(iss, line);
    // Line 3: NSS, NSX_
    std::getline(iss, line);
    int NSX_local = 0;
    {
        std::istringstream ls(line);
        ls >> state.NSS >> NSX_local;
    }
    // Line 4: column headers
    std::getline(iss, line);

    for (int K = 0; K < state.NSS; ++K) {
        // Read type header: title, RSS, GSS, DSS, TSS, WSS
        std::getline(iss, line);
        state.TITLSS[K] = line.substr(0, 10);
        std::istringstream ls(line.substr(10));
        ls >> state.RSS[K] >> state.GSS[K] >> state.DSS[K]
           >> state.TSS[K] >> state.WSS[K];

        // Read wavelength data from bin 5 to NSX_ (Fortran J=5,NSX_)
        for (int J = 4; J < NSX_local; ++J) {
            std::getline(iss, line);
            std::istringstream wls(line);
            int JSS;
            double WJSS, XNDR, XNDI;
            wls >> JSS >> WJSS >> XNDR >> XNDI;
            wls >> state.QSS[J][K] >> state.SSS[J][K];
            for (int I = 1; I < 8; ++I) {
                wls >> state.PSS[I][J][K];
            }
            state.PSS[0][J][K] = 1.0;
        }
    }

    // Replicate SSA data for J=0:3 from J=4
    for (int K = 0; K < state.NSS; ++K) {
        for (int J = 0; J < 4; ++J) {
            state.QSS[J][K] = state.QSS[4][K];
            state.SSS[J][K] = state.SSS[4][K];
            for (int I = 0; I < 8; ++I) {
                state.PSS[I][J][K] = state.PSS[I][4][K];
            }
        }
    }
}

// =========================================================================
// RD_MIE: Read aerosol Mie scattering data
// =========================================================================

inline void RD_MIE(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_MIE (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "FJX_scat-aer.dat", Tables::FJX_scat_aer));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2-3: column header lines
    std::getline(iss, line);
    std::getline(iss, line);

    for (int J = 0; J < A_; ++J) {
        std::getline(iss, line);
        if (line.empty() || iss.eof()) break;

        // Format: "  04|1x,a12,1x,2f7.3,1x,a120"
        // Parse JAA from first 4 chars
        std::istringstream hdr(line);
        int JAA;
        // Try to read the integer (before the '|')
        std::string jaa_str = line.substr(0, 4);
        // Trim whitespace
        size_t start = jaa_str.find_first_not_of(' ');
        if (start == std::string::npos) break;
        JAA = std::stoi(jaa_str.substr(start));

        if (JAA <= 0) break;

        // Title is chars 5..16 (after '|')
        state.TITLAA[J] = line.substr(5, 12);
        // RAA and DAA
        std::istringstream vals(line.substr(18, 14));
        vals >> state.RAA[J] >> state.DAA[J];

        // Read 5 wavelength lines
        for (int K = 0; K < 5; ++K) {
            std::getline(iss, line);
            std::istringstream wls(line);
            double waa;
            wls >> waa >> state.QAA[K][J] >> state.SAA[K][J];
            state.WAA[K][J] = waa;
            for (int I = 1; I < 8; ++I) {
                wls >> state.PAA[I][K][J];
            }
            state.PAA[0][K][J] = 1.0;
        }
        state.NAA = J + 1;
    }
}

// =========================================================================
// RD_UM: Read UMichigan aerosol optical data
// =========================================================================

inline void RD_UM(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_UM (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "FJX_scat-UMa.dat", Tables::FJX_scat_UMa));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2: wavelengths
    std::getline(iss, line);
    {
        std::istringstream ls(line.substr(5));  // skip "wavel" prefix
        for (int i = 0; i < 6; ++i) {
            ls >> state.WMM[i];
        }
    }

    // 33 aerosol types
    for (int L = 0; L < 33; ++L) {
        // Type title line (4 chars)
        std::getline(iss, line);
        // 21 relative humidity values
        for (int K = 0; K < 21; ++K) {
            std::getline(iss, line);
            std::istringstream ls(line);
            // 6 wavelengths x 3 optic vars = 18 values per line
            for (int J = 0; J < 6; ++J) {
                for (int I = 0; I < 3; ++I) {
                    ls >> state.UMAER[I][J][K][L];
                }
            }
        }
    }

    // Collapse wavelengths: drop 550nm (index 3), shift 600 and 1000
    state.WMM[3] = state.WMM[4];
    state.WMM[4] = state.WMM[5];
    for (int L = 0; L < 33; ++L) {
        for (int K = 0; K < 21; ++K) {
            for (int I = 0; I < 3; ++I) {
                state.UMAER[I][3][K][L] = state.UMAER[I][4][K][L];
                state.UMAER[I][4][K][L] = state.UMAER[I][5][K][L];
            }
        }
    }
}

// =========================================================================
// RD_GEO: Read GEOMIP SSA scattering data
// =========================================================================

inline void RD_GEO(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_GEO (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "FJX_scat-geo.dat", Tables::FJX_scat_geo));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2: NGG
    std::getline(iss, line);
    state.NGG = std::stoi(line);
    // Line 3-4: column headers
    std::getline(iss, line);
    std::getline(iss, line);

    for (int K = 0; K < state.NGG; ++K) {
        // Read type header with R, G, dens, T, wt
        std::getline(iss, line);
        std::istringstream hdr(line.substr(10));
        double G1, G2, G3;
        hdr >> state.RGG[K] >> G1 >> state.DGG[K] >> G2 >> G3;

        // Read wavelength data from bin 5 to 27 (Fortran J=5,27)
        for (int J = 4; J < 27; ++J) {
            std::getline(iss, line);
            std::istringstream wls(line.substr(2));
            double WGGJ, XNDR, XNDI;
            wls >> WGGJ >> XNDR >> XNDI;
            wls >> state.QGG[J][K] >> state.SGG[J][K];
            for (int I = 1; I < 8; ++I) {
                wls >> state.PGG[I][J][K];
            }
            state.PGG[0][J][K] = 1.0;
        }
    }

    // Replicate for J=0:3 from J=4
    for (int K = 0; K < state.NGG; ++K) {
        for (int J = 0; J < 4; ++J) {
            state.QGG[J][K] = state.QGG[4][K];
            state.SGG[J][K] = state.SGG[4][K];
            for (int I = 0; I < 8; ++I) {
                state.PGG[I][J][K] = state.PGG[I][4][K];
            }
        }
    }
}

// =========================================================================
// RD_PROF: Read T and O3 reference profiles from atmos_std table
// =========================================================================

inline void RD_PROF(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_PROF (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "atmos_std.dat", Tables::atmos_std));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2: NTLATS, NTMONS
    std::getline(iss, line);
    int NTLATS = 0, NTMONS = 0;
    {
        std::istringstream ls(line);
        ls >> NTLATS >> NTMONS;
    }

    int N216 = std::min(216, NTLATS * NTMONS);
    for (int IA = 0; IA < N216; ++IA) {
        // Read LAT, MON header
        std::getline(iss, line);
        int LAT = 0, MON = 0;
        // Format: "A-85LAT 1M" -> parse latitude from position and month
        // Actually format is "1X,I3,3X,I2" in Fortran
        std::istringstream ls(line.substr(1));
        ls >> LAT;
        std::string rest = line.substr(4);
        std::istringstream ls2(rest.substr(3));
        ls2 >> MON;

        int M = std::min(11, std::max(0, MON - 1));  // 0-indexed month
        int L = std::min(17, std::max(0, (LAT + 95) / 10 - 1));  // 0-indexed latitude

        // Read T_REF: 41 values over multiple lines (11 per line)
        {
            int read = 0;
            while (read < 41) {
                std::getline(iss, line);
                std::istringstream tls(line.substr(3));  // skip "T1 " prefix
                double val;
                while (tls >> val && read < 41) {
                    state.T_REF[read][L][M] = val;
                    read++;
                }
            }
        }

        // Read O_REF: 31 values
        {
            int read = 0;
            while (read < 31) {
                std::getline(iss, line);
                std::istringstream ols(line.substr(3));  // skip "Z1 " prefix
                double val;
                while (ols >> val && read < 31) {
                    state.O_REF[read][L][M] = val;
                    read++;
                }
            }
        }
    }

    // Extend climatology to 100 km (LREF=51 layers)
    double OFAC = std::exp(-2.0e5 / 5.0e5);
    for (int I = 31; I < LREF; ++I) {
        double OFAK = std::pow(OFAC, I - 30);
        for (int M = 0; M < NTMONS; ++M) {
            for (int L2 = 0; L2 < NTLATS; ++L2) {
                state.O_REF[I][L2][M] = state.O_REF[30][L2][M] * OFAK;
            }
        }
    }
    for (int L2 = 0; L2 < NTLATS; ++L2) {
        for (int M = 0; M < NTMONS; ++M) {
            for (int I = 41; I < LREF; ++I) {
                state.T_REF[I][L2][M] = state.T_REF[40][L2][M];
            }
        }
    }
}

// =========================================================================
// RD_TRPROF: Read H2O and CH4 reference profiles from atmos_h2och4 table
// =========================================================================

inline void RD_TRPROF(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_TRPROF (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "atmos_h2och4.dat", Tables::atmos_h2och4));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2: NTLATS, NTMONS
    std::getline(iss, line);
    int NTLATS = 0, NTMONS = 0;
    {
        std::istringstream ls(line);
        ls >> NTLATS >> NTMONS;
    }

    int N216 = std::min(216, NTLATS * NTMONS);
    for (int IA = 0; IA < N216; ++IA) {
        // Read LAT, MON header
        std::getline(iss, line);
        int LAT = 0, MON = 0;
        std::istringstream ls(line.substr(1));
        ls >> LAT;
        std::string rest = line.substr(4);
        std::istringstream ls2(rest.substr(3));
        ls2 >> MON;

        int M = std::min(11, std::max(0, MON - 1));
        int L = std::min(17, std::max(0, (LAT + 95) / 10 - 1));

        // Read H2O_REF: 31 values over multiple lines
        {
            int read = 0;
            while (read < 31) {
                std::getline(iss, line);
                std::istringstream hls(line.substr(3));
                double val;
                while (hls >> val && read < 31) {
                    state.H2O_REF[read][L][M] = val;
                    read++;
                }
            }
        }

        // Read CH4_REF: 31 values over multiple lines
        {
            int read = 0;
            while (read < 31) {
                std::getline(iss, line);
                std::istringstream cls(line.substr(3));
                double val;
                while (cls >> val && read < 31) {
                    state.CH4_REF[read][L][M] = val;
                    read++;
                }
            }
        }
    }

    // Extend climatology to LREF
    for (int L2 = 0; L2 < NTLATS; ++L2) {
        for (int M = 0; M < NTMONS; ++M) {
            for (int I = 31; I < LREF; ++I) {
                state.H2O_REF[I][L2][M] = state.H2O_REF[30][L2][M];
                state.CH4_REF[I][L2][M] = state.CH4_REF[30][L2][M];
            }
        }
    }
}

// =========================================================================
// RD_JS_JX: Read photolysis rate mapping data from FJX_j2j table
// =========================================================================

inline void RD_JS_JX(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_JS_JX (C++ port)";
    rc = CLDJ_SUCCESS;

    // Initialize mapping arrays
    for (int i = 0; i < JVN_; ++i) {
        state.JLABEL[i] = "------";
        state.JVMAP[i] = "------";
        state.JFACTA[i] = 0.0;
    }

    std::istringstream iss(load_table_or_default(datadir, "FJX_j2j.dat", Tables::FJX_j2j));
    std::string line;

    // Line 1: header/title
    std::getline(iss, line);

    // Read entries until JJ == 9999
    for (int entry = 0; entry < JVN_; ++entry) {
        std::getline(iss, line);
        if (line.empty() || iss.eof()) break;

        // Format: "i4,1x,a50,4x,f5.3,2x,a6"
        // Example: "   1 NO        PHOTON    N         O                       1.000 /NO    /"
        int JJ = 0;
        try {
            JJ = std::stoi(line.substr(0, 4));
        } catch (...) {
            break;
        }

        if (JJ == 9999) break;

        if (JJ > JVN_) {
            CLOUDJ_ERROR("JVN_ must be >= number of entries in FJX_j2j", thisloc, rc);
            return;
        }

        // 0-indexed JJ
        int idx = JJ - 1;

        // Reaction label (50 chars starting at position 5)
        std::string T_REACT = line.substr(5, 50);
        state.JLABEL[idx] = T_REACT;

        // Factor (5 chars at position 59)
        double F_FJX = 0.0;
        try {
            F_FJX = std::stod(line.substr(59, 5));
        } catch (...) {
            F_FJX = 1.0;
        }
        state.JFACTA[idx] = F_FJX;

        // FJX name (6 chars at position 66, between '/' delimiters)
        std::string T_FJX = "------";
        size_t slash1 = line.find('/', 55);
        if (slash1 != std::string::npos && slash1 + 7 <= line.size()) {
            T_FJX = line.substr(slash1 + 1, 6);
        }
        state.JVMAP[idx] = T_FJX;
        state.NRATJ = JJ;

        // Extract reaction name (first 10 chars of T_REACT)
        state.RNAMES[idx] = T_REACT.substr(0, 10);
        // Trim trailing spaces
        size_t end = state.RNAMES[idx].find_last_not_of(' ');
        if (end != std::string::npos)
            state.RNAMES[idx] = state.RNAMES[idx].substr(0, end + 1);

        // Compute branch number
        state.BRANCH[idx] = 1;
        for (int K = 0; K < idx; ++K) {
            if (state.RNAMES[idx] == state.RNAMES[K]) {
                state.BRANCH[idx] = state.BRANCH[K] + 1;
            }
        }
    }

    // Map JVMAP onto TITLEJX to set JIND
    for (int K = 0; K < state.NRATJ; ++K) {
        state.JIND[K] = 0;
        for (int J = 0; J < state.NJX; ++J) {
            if (state.JVMAP[K] == state.TITLEJX[J]) {
                state.JIND[K] = J + 1;  // 1-based index for compatibility
            }
        }
    }
}

// =========================================================================
// RD_SSAPROF: Read SSA-GEO reference profiles from atmos_geomip table
// =========================================================================

inline void RD_SSAPROF(const std::string& datadir, CloudJState& state, int& rc) {
    const std::string thisloc = " -> at RD_SSAPROF (C++ port)";
    rc = CLDJ_SUCCESS;

    std::istringstream iss(load_table_or_default(datadir, "atmos_geomip.dat", Tables::atmos_geomip));
    std::string line;

    // Line 1: title
    std::getline(iss, line);
    // Line 2: comment
    std::getline(iss, line);
    // Line 3: "Pressure levels 11:28" header
    std::getline(iss, line);
    // Line 4: P_GREF values (19 values)
    std::getline(iss, line);
    {
        std::istringstream ls(line);
        for (int L = 0; L < 19; ++L) {
            ls >> state.P_GREF[L];
        }
    }
    // Line 5: "Pressure Mass" header
    std::getline(iss, line);
    // Line 6: pressure mass values (skip)
    std::getline(iss, line);
    // Line 7: "Latitudes 1:32" header
    std::getline(iss, line);
    // Line 8: Y_GREF[0:31] (32 values)
    std::getline(iss, line);
    {
        std::istringstream ls(line);
        for (int L = 0; L < 32; ++L) {
            ls >> state.Y_GREF[L];
        }
    }
    // Line 9: blank/header
    std::getline(iss, line);
    // Line 10: Y_GREF[63:32] (reversed, 32 values)
    std::getline(iss, line);
    {
        std::istringstream ls(line);
        for (int L = 63; L >= 32; --L) {
            ls >> state.Y_GREF[L];
        }
    }

    // Read R_GREF section
    std::getline(iss, line);  // "R = R-effective..." header
    for (int M = 0; M < 12; ++M) {
        std::getline(iss, line);  // month header line
        for (int J = 0; J < 64; ++J) {
            std::getline(iss, line);
            // Skip first 11 chars (index info), then read 18 values
            std::istringstream ls(line.substr(11));
            for (int L = 0; L < 18; ++L) {
                ls >> state.R_GREF[J][L][M];
            }
            state.R_GREF[J][18][M] = 0.0;
        }
    }

    // Read X_GREF section
    std::getline(iss, line);  // "X = mass fraction..." header
    for (int M = 0; M < 12; ++M) {
        std::getline(iss, line);  // month header line
        for (int J = 0; J < 64; ++J) {
            std::getline(iss, line);
            std::istringstream ls(line.substr(11));
            for (int L = 0; L < 18; ++L) {
                ls >> state.X_GREF[J][L][M];
            }
            state.X_GREF[J][18][M] = 0.0;
        }
    }

    // Read A_GREF section
    std::getline(iss, line);  // "A = 4*pi*R^2..." header
    for (int M = 0; M < 12; ++M) {
        std::getline(iss, line);  // month header line
        for (int J = 0; J < 64; ++J) {
            std::getline(iss, line);
            std::istringstream ls(line.substr(11));
            for (int L = 0; L < 18; ++L) {
                ls >> state.A_GREF[J][L][M];
            }
            state.A_GREF[J][18][M] = 0.0;
        }
    }
}

// =========================================================================
// INIT_CLDJ: Main initialization driver for Cloud-J
// Ported from cldj_init_mod.F90 INIT_CLDJ subroutine
// =========================================================================

inline void INIT_CLDJ(
    bool AMIROOT,
    const std::string& DATADIR,
    int NLEVELS,
    int NLEVELS_WITH_CLOUD,
    std::vector<std::string>& TITLEJXX,
    int NJXU,
    double ATAU_in,
    double ATAU0_in,
    int NWBIN_in,
    int CLDFLAG_in,
    double CLDCOR_in,
    int LNRG_in,
    int ATM0_in,
    bool use_H2O_UV_abs,
    int& NJXX,
    CloudJState& state,
    int& rc)
{
    const std::string thisloc = " -> at INIT_CLDJ (C++ port)";
    rc = CLDJ_SUCCESS;

    // ---------------------------------------------------------------
    // Set grid dimension parameters (non-standalone only)
    // ---------------------------------------------------------------
#ifndef MODEL_STANDALONE
    state.L_     = NLEVELS;
    state.L1_    = NLEVELS + 1;
    state.L2_    = NLEVELS + 2;
    state.LWEPAR = NLEVELS_WITH_CLOUD;
#else
    // Suppress unused parameter warnings in standalone mode
    (void)NLEVELS;
    (void)NLEVELS_WITH_CLOUD;
#endif

    // ---------------------------------------------------------------
    // Set runtime configuration from input parameters
    // ---------------------------------------------------------------
    state.ATAU     = ATAU_in;
    state.ATAU0    = ATAU0_in;
    state.CLDCOR   = CLDCOR_in;
    state.NWBIN    = NWBIN_in;
    state.LNRG     = LNRG_in;
    state.ATM0     = ATM0_in;
    state.CLDFLAG  = CLDFLAG_in;
    state.USEH2OUV = use_H2O_UV_abs;

    // v7.7 safety fix: if LNRG != 6, force CLDCOR = 0
    if (state.LNRG != 6) {
        state.CLDCOR = 0.0;
    }

    // ---------------------------------------------------------------
    // Validate W_ == 18
    // ---------------------------------------------------------------
    if (W_ != 18) {
        CLOUDJ_ERROR("Invalid no. wavelengths: W_ != 18", thisloc, rc);
        return;
    }

    // ---------------------------------------------------------------
    // Initialize NSJSUB and SJSUB defaults (Cloud-J: no sub-bins)
    // ---------------------------------------------------------------
    for (int I = 0; I < S_; ++I) {
        state.NSJSUB[I] = NGC[I];
    }
    for (int I = 0; I < SX_; ++I) {
        state.SJSUB[I][0] = 1.0;
        for (int J = 1; J < 16; ++J) {
            state.SJSUB[I][J] = 0.0;
        }
    }

    // ---------------------------------------------------------------
    // Call RD_XXX: read spectral data
    // ---------------------------------------------------------------
    RD_XXX(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_XXX", thisloc, rc);
        return;
    }

    // Set up KDOKR / LDOKR index arrays
    for (int i = 0; i < 100; ++i) {
        state.KDOKR[i] = 0;
        state.LDOKR[i] = 0;
    }
    int KR = 0;
    for (int K = 0; K < S_; ++K) {
        for (int J = 0; J < state.NSJSUB[K]; ++J) {
            if (KR >= 100) break;
            state.KDOKR[KR] = K + 1;  // 1-based for Fortran compat
            if (state.FL[K] > 0.0) {
                state.LDOKR[KR] = 1;
            } else {
                state.LDOKR[KR] = 0;
            }
            KR++;
        }
    }
    if (KR != W_ + W_r) {
        CLOUDJ_ERROR("Error with sub-bin setup: KDOKR", thisloc, rc);
        return;
    }

    // ---------------------------------------------------------------
    // Call RD_CLD: read cloud scattering data
    // ---------------------------------------------------------------
    RD_CLD(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_CLD", thisloc, rc);
        return;
    }

    // ---------------------------------------------------------------
    // Call RD_SSA: read stratospheric sulfate aerosol data
    // ---------------------------------------------------------------
    RD_SSA(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_SSA", thisloc, rc);
        return;
    }

    // ---------------------------------------------------------------
    // Call RD_MIE: read aerosol Mie scattering data
    // ---------------------------------------------------------------
    RD_MIE(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_MIE", thisloc, rc);
        return;
    }

    // ---------------------------------------------------------------
    // Call RD_UM: read UMich aerosol data (not used in GEOS-Chem)
    // ---------------------------------------------------------------
#ifdef MODEL_GEOSCHEM
    for (int i = 0; i < 6; ++i) state.WMM[i] = 0.0;
    std::fill(&state.UMAER[0][0][0][0],
              &state.UMAER[0][0][0][0] + sizeof(state.UMAER)/sizeof(double), 0.0);
#else
    RD_UM(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_UM", thisloc, rc);
        return;
    }
#endif

    // ---------------------------------------------------------------
    // Call RD_GEO: read GEOMIP aerosol scattering data
    // ---------------------------------------------------------------
#ifdef MODEL_GEOSCHEM
    state.NGG = 0;
    std::fill(&state.RGG[0], &state.RGG[0] + GGA_, 0.0);
    std::fill(&state.DGG[0], &state.DGG[0] + GGA_, 0.0);
    std::fill(&state.QGG[0][0], &state.QGG[0][0] + SX_ * GGA_, 0.0);
    std::fill(&state.SGG[0][0], &state.SGG[0][0] + SX_ * GGA_, 0.0);
    std::fill(&state.PGG[0][0][0], &state.PGG[0][0][0] + 8 * SX_ * GGA_, 0.0);
#else
    RD_GEO(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_GEO", thisloc, rc);
        return;
    }
#endif

    // ---------------------------------------------------------------
    // Call RD_PROF / RD_TRPROF: read climatology profiles
    // ---------------------------------------------------------------
#ifdef MODEL_STANDALONE
    RD_PROF(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_PROF", thisloc, rc);
        return;
    }

    RD_TRPROF(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_TRPROF", thisloc, rc);
        return;
    }
#else
    // Zero out climatology arrays when not in standalone mode
    std::fill(&state.T_REF[0][0][0],
              &state.T_REF[0][0][0] + LREF * JREF * 12, 0.0);
    std::fill(&state.O_REF[0][0][0],
              &state.O_REF[0][0][0] + LREF * JREF * 12, 0.0);
    std::fill(&state.H2O_REF[0][0][0],
              &state.H2O_REF[0][0][0] + LREF * JREF * 12, 0.0);
    std::fill(&state.CH4_REF[0][0][0],
              &state.CH4_REF[0][0][0] + LREF * JREF * 12, 0.0);
#endif

    // ---------------------------------------------------------------
    // Call RD_SSAPROF: read GeoMIP SSA reference profiles
    // ---------------------------------------------------------------
#ifdef MODEL_GEOSCHEM
    std::fill(&state.R_GREF[0][0][0],
              &state.R_GREF[0][0][0] + 64 * LGREF * 12, 0.0);
    std::fill(&state.X_GREF[0][0][0],
              &state.X_GREF[0][0][0] + 64 * LGREF * 12, 0.0);
    std::fill(&state.A_GREF[0][0][0],
              &state.A_GREF[0][0][0] + 64 * LGREF * 12, 0.0);
    std::fill(&state.Y_GREF[0], &state.Y_GREF[0] + 64, 0.0);
    std::fill(&state.P_GREF[0], &state.P_GREF[0] + LGREF, 0.0);
#else
    RD_SSAPROF(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_SSAPROF", thisloc, rc);
        return;
    }
#endif

    // ---------------------------------------------------------------
    // Populate TITLEJXX output and set NJXX from NJX
    // ---------------------------------------------------------------
    NJXX = state.NJX;
    TITLEJXX.resize(std::max(NJXU, state.NJX));
    for (int J = 0; J < state.NJX; ++J) {
        TITLEJXX[J] = state.TITLEJX[J];
    }

    // ---------------------------------------------------------------
    // Call RD_JS_JX: read photolysis rate mapping
    // ---------------------------------------------------------------
    RD_JS_JX(DATADIR, state, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RD_JS_JX", thisloc, rc);
        return;
    }

    // ---------------------------------------------------------------
    // Call RANSET: generate pseudo-random number sequence
    // ---------------------------------------------------------------
    int RANSEED = 66;
    RANSET(state.RAN4, RANSEED, rc);
    if (rc != CLDJ_SUCCESS) {
        CLOUDJ_ERROR("Error in RANSET", thisloc, rc);
        return;
    }
}

} // namespace CloudJ::Init

#endif // CLOUDJ_INIT_HPP
