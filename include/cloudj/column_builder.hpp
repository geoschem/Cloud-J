#ifndef CLOUDJ_COLUMN_BUILDER_HPP
#define CLOUDJ_COLUMN_BUILDER_HPP

// Shared atmospheric-column derivation.
//
// The standalone reference driver reads a climatological column description
// (tables/atmos_PTClds.dat), turns it into the full set of arrays that
// CLOUD_JX consumes (pressure/height edges, temperature, air/O3/CH4/H2O
// columns, relative humidity, aerosol loading, cloud fraction/water paths and
// effective radii, and the surface-reflectivity table), and then calls the
// solver. This header factors that one-to-one derivation out of the driver so
// that the high-level convenience API can build exactly the same column and
// obtain Fortran-identical J-values for the same physical input, rather than
// maintaining a second, independently-written copy of the formulas.
//
// The routines here mirror the standalone driver operation-for-operation; any
// difference in the resulting J-values would be a bug in one caller, not a
// divergence between two implementations.

#include <algorithm>
#include <cloudj/osa.hpp>
#include <cloudj/photo_jx.hpp>
#include <cloudj/state.hpp>
#include <cmath>
#include <cstring>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

namespace CloudJ {
namespace Column {

// Raw column description as read from atmos_PTClds.dat. These are the inputs a
// host model would supply; everything derived from them lives in Derived.
struct Inputs {
  int month = 0;
  int ilat = 0;
  double psurf = 0.0;
  double albedo[5] = {};
  double wind = 0.0;
  double chlr = 0.0;

  // Per-level tables. Edge coefficients are sized L2_ so the top two entries
  // (index L1_ and L2_-1) can be zero-filled the way the reference does; the
  // remaining per-level arrays are sized L1_.
  std::vector<double> etaa;
  std::vector<double> etab;
  std::vector<double> tinp;
  std::vector<double> rhinp;
  std::vector<double> aer1;
  std::vector<double> aer2;
  std::vector<int> naa1;
  std::vector<int> naa2;

  // Per-cloud-level tables, sized LWEPAR.
  std::vector<double> cldfrw;
  std::vector<double> cldlwcw;
  std::vector<double> cldiwcw;

  Inputs() {
    constexpr int L1 = CloudJState::L1_;
    constexpr int L2 = CloudJState::L2_;
    constexpr int LWE = CloudJState::LWEPAR;
    etaa.assign(L2, 0.0);
    etab.assign(L2, 0.0);
    tinp.assign(L1, 0.0);
    rhinp.assign(L1, 0.0);
    aer1.assign(L1, 0.0);
    aer2.assign(L1, 0.0);
    naa1.assign(L1, 0);
    naa2.assign(L1, 0);
    cldfrw.assign(LWE, 0.0);
    cldlwcw.assign(LWE, 0.0);
    cldiwcw.assign(LWE, 0.0);
  }
};

// Fully derived column: every array CLOUD_JX needs, ready to be handed to the
// solver for a given solar geometry.
struct Derived {
  std::vector<double> ppp;   // edge pressure (hPa), size L2_
  std::vector<double> zzz;   // edge altitude (cm), size L2_
  std::vector<double> ttt;   // temperature (K), size L1_
  std::vector<double> ddd;   // air column density, size L1_
  std::vector<double> ooo;   // O3 column, size L1_
  std::vector<double> ccc;   // CH4 column, size L1_
  std::vector<double> hhh;   // H2O column, size L1_
  std::vector<double> rrr;   // relative humidity, size L1_
  std::vector<double> aersp; // aerosol optical path, size L1_ * AN_
  std::vector<int> ndxaer;   // aerosol type index, size L1_ * AN_
  std::vector<double> clf;   // cloud fraction (mutated by CLOUD_JX), size L1_
  std::vector<int> cldiw;    // cloud phase index, size L1_
  std::vector<double> lwp;   // liquid water path, size L1_
  std::vector<double> iwp;   // ice water path, size L1_
  std::vector<double> reffl; // liquid effective radius, size L1_
  std::vector<double> reffi; // ice effective radius, size L1_
  std::vector<double> clf0;  // pristine cloud fraction to reset clf each call
  int ltop = 0;              // number of levels that can carry cloud

  Derived() {
    constexpr int L1 = CloudJState::L1_;
    constexpr int L2 = CloudJState::L2_;
    ppp.assign(L2, 0.0);
    zzz.assign(L2, 0.0);
    ttt.assign(L1, 0.0);
    ddd.assign(L1, 0.0);
    ooo.assign(L1, 0.0);
    ccc.assign(L1, 0.0);
    hhh.assign(L1, 0.0);
    rrr.assign(L1, 0.0);
    aersp.assign(L1 * AN_, 0.0);
    ndxaer.assign(L1 * AN_, 0);
    clf.assign(L1, 0.0);
    cldiw.assign(L1, 0);
    lwp.assign(L1, 0.0);
    iwp.assign(L1, 0.0);
    reffl.assign(L1, 0.0);
    reffi.assign(L1, 0.0);
    clf0.assign(L1, 0.0);
    ltop = CloudJState::LWEPAR;
  }
};

// Parse the atmos_PTClds.dat fixed-format profile from an already-open stream
// (used by tests that embed the reference column without touching the
// filesystem).
inline void parse_stream(std::istream& infile, Inputs& out) {
  constexpr int L1 = CloudJState::L1_;
  std::string line;

  // Line 1: title (skip).
  std::getline(infile, line);
  // Line 2: MONTH, ILAT (format 2i5).
  std::getline(infile, line);
  out.month = std::stoi(line.substr(0, 5));
  out.ilat = std::stoi(line.substr(5, 5));
  // Line 3: PSURF (format f5.0).
  std::getline(infile, line);
  out.psurf = std::stod(line.substr(0, 5));
  // Line 4: ALBEDO(5) for the incident ray (format f5.2).
  std::getline(infile, line);
  out.albedo[4] = std::stod(line.substr(0, 5));
  // Line 5: ALBEDO(1:4) for the four quadrature angles (format 4f5.2).
  std::getline(infile, line);
  for (int i = 0; i < 4; ++i) {
    out.albedo[i] = std::stod(line.substr(i * 5, 5));
  }
  // Line 6: surface wind and chlorophyll for the ocean-surface albedo.
  std::getline(infile, line);
  {
    std::istringstream ss(line);
    ss >> out.wind >> out.chlr;
  }
  // Line 7: column header (skip).
  std::getline(infile, line);

  // L1 levels of per-level atmosphere data. ZOFL is parsed for completeness but
  // is not consumed downstream, matching the reference.
  for (int L = 0; L < L1; ++L) {
    std::getline(infile, line);
    std::istringstream ss(line);
    int idx;
    double zofl;
    ss >> idx >> out.etaa[L] >> out.etab[L] >> out.tinp[L] >> out.rhinp[L]
       >> zofl >> out.aer1[L] >> out.naa1[L] >> out.aer2[L] >> out.naa2[L];
  }

  // Cloud table header (skip).
  std::getline(infile, line);

  // LWEPAR cloud levels, stored reversed (the file lists them top-down).
  for (int L = CloudJState::LWEPAR - 1; L >= 0; --L) {
    std::getline(infile, line);
    std::istringstream ss(line);
    int idx;
    ss >> idx >> out.cldfrw[L] >> out.cldlwcw[L] >> out.cldiwcw[L];
  }
}

// Parse the atmos_PTClds.dat fixed-format profile from a file path, matching
// the reference reader field-for-field. Returns false if the file cannot be
// opened.
inline bool parse_atmos_ptclds(const std::string& path, Inputs& out) {
  std::ifstream infile(path);
  if (!infile.is_open()) {
    return false;
  }
  parse_stream(infile, out);
  return true;
}

// Build the full derived column from parsed inputs. The supplied state is used
// for the Fast-JX climatology (O3, temperature, CH4) exactly as the reference
// driver does. This reproduces the standalone derivation step for step.
inline void build_column(const Inputs& in, const CloudJState& state,
                         Derived& out) {
  constexpr int L1 = CloudJState::L1_;
  constexpr int L2 = CloudJState::L2_;
  constexpr int LWE = CloudJState::LWEPAR;

  // Edge pressure: PPP = eta-A + eta-B * PSURF, with the two top entries zero.
  std::vector<double> etaa = in.etaa;
  std::vector<double> etab = in.etab;
  etaa[L2 - 1] = 0.0;
  etab[L2 - 1] = 0.0;
  for (int L = 0; L < L2; ++L) {
    out.ppp[L] = etaa[L] + etab[L] * in.psurf;
  }

  // Fast-JX climatology for this month/latitude on the pressure grid.
  std::vector<double> o3mix(L1, 0.0);
  std::vector<double> ch4mix(L1, 0.0);
  double ylat = static_cast<double>(in.ilat);
  PhotoJX::ACLIM_FJX(in.month, ylat, out.ppp.data(), out.ttt.data(),
                     o3mix.data(), ch4mix.data(), L1, state);

  // The file's own temperature and relative humidity override the climatology;
  // O3 and CH4 keep the climatological values.
  out.rrr = in.rhinp;
  for (int L = 0; L < L1; ++L) {
    out.ttt[L] = in.tinp[L];
  }

  // Altitudes, air density, and O3/CH4 columns from the pressure grid.
  out.zzz[0] = 16.0e5 * std::log10(1013.25 / out.ppp[0]);
  for (int L = 0; L < CloudJState::L_; ++L) {
    out.ddd[L] = (out.ppp[L] - out.ppp[L + 1]) * MASFAC;
    double scaleh = 1.3806e-19 * MASFAC * out.ttt[L];
    out.zzz[L + 1] = out.zzz[L] - (std::log(out.ppp[L + 1] / out.ppp[L]) * scaleh);
    out.ooo[L] = out.ddd[L] * o3mix[L] * 1.0e-6;
    out.ccc[L] = out.ddd[L] * ch4mix[L] * 1.0e-9;
  }
  {
    // Extra top layer above the CTM: fixed scale-height offset, no temperature
    // integration.
    int L = CloudJState::L_;
    out.zzz[L + 1] = out.zzz[L] + ZZHT;
    out.ddd[L] = (out.ppp[L] - out.ppp[L + 1]) * MASFAC;
    out.ooo[L] = out.ddd[L] * o3mix[L] * 1.0e-6;
    out.ccc[L] = out.ddd[L] * ch4mix[L] * 1.0e-9;
  }

  // H2O column: exponentially decaying mixing ratio with a fixed floor.
  for (int L = 0; L < L1; ++L) {
    out.hhh[L] = out.ddd[L] * std::max(0.030 * std::exp(-out.zzz[L] / 2.2e5),
                                       2.0e-6);
  }

  // Aerosol loading: two species per level, column-major [L][AN_].
  for (int L = 0; L < CloudJState::L_; ++L) {
    out.ndxaer[L + L1 * 0] = in.naa1[L];
    out.aersp[L + L1 * 0] = in.aer1[L];
    out.ndxaer[L + L1 * 1] = in.naa2[L];
    out.aersp[L + L1 * 1] = in.aer2[L];
  }

  // Cloud processing. Detect an entirely clear column first, then build the
  // per-level fraction, phase, water paths, and effective radii.
  out.ltop = LWE;
  double max_clf = 0.0;
  for (int L = 0; L < LWE; ++L) {
    max_clf = std::max(max_clf, in.cldfrw[L]);
  }
  if (max_clf <= 0.005) {
    std::memset(out.iwp.data(), 0, L1 * sizeof(double));
    std::memset(out.reffi.data(), 0, L1 * sizeof(double));
    std::memset(out.lwp.data(), 0, L1 * sizeof(double));
    std::memset(out.reffl.data(), 0, L1 * sizeof(double));
  }

  std::vector<double> wlc(L1, 0.0);
  std::vector<double> wic(L1, 0.0);
  for (int L = 0; L < out.ltop; ++L) {
    out.cldiw[L] = 0;
    double cf = in.cldfrw[L];
    if (cf > 0.005) {
      out.clf[L] = cf;
      wlc[L] = in.cldlwcw[L] / cf;
      wic[L] = in.cldiwcw[L] / cf;
      if (wlc[L] > 1.0e-11) out.cldiw[L] = 1;
      if (wic[L] > 1.0e-11) out.cldiw[L] = out.cldiw[L] + 2;
    } else {
      out.clf[L] = 0.0;
      wlc[L] = 0.0;
      wic[L] = 0.0;
    }
  }

  for (int L = 0; L < out.ltop; ++L) {
    if (wic[L] > 1.0e-12) {
      double pdel = out.ppp[L] - out.ppp[L + 1];
      double zdel = (out.zzz[L + 1] - out.zzz[L]) * 0.01; // m
      out.iwp[L] = 1000.0 * wic[L] * pdel * G100;         // g/m2
      double icwc = out.iwp[L] / zdel;                    // g/m3
      out.reffi[L] = 164.0 * std::pow(icwc, 0.23);
    } else {
      out.iwp[L] = 0.0;
      out.reffi[L] = 0.0;
    }
    if (wlc[L] > 1.0e-12) {
      double pmid = 0.5 * (out.ppp[L] + out.ppp[L + 1]);
      double pdel = out.ppp[L] - out.ppp[L + 1];
      double f1 = 0.005 * (pmid - 610.0);
      f1 = std::min(1.0, std::max(0.0, f1));
      out.lwp[L] = 1000.0 * wlc[L] * pdel * G100; // g/m2
      out.reffl[L] = 9.6 * f1 + 12.68 * (1.0 - f1);
    } else {
      out.lwp[L] = 0.0;
      out.reffl[L] = 0.0;
    }
  }

  // Keep a pristine copy of the cloud fraction: CLOUD_JX consumes CLDF in
  // place, so callers reset from this before each solve.
  out.clf0 = out.clf;
}

// Reset the (solver-mutated) cloud fraction back to its pristine values. Call
// this before every CLOUD_JX invocation that reuses a Derived column.
inline void reset_cloud_fraction(Derived& col) {
  for (int L = 0; L < col.ltop; ++L) {
    col.clf[L] = col.clf0[L];
  }
}

// Build the surface-reflectivity table for a given cos(zenith). The reference
// computes the ocean-surface albedo and then overrides every entry with the
// read-in broadband albedo, so the result is the albedo broadcast across the
// spectral bins in the [angle + 5*bin] layout CLOUD_JX expects.
inline void build_rfl(const CloudJState& state, double u0,
                      const double albedo[5], double wind, double chlr,
                      std::vector<double>& rfl) {
  constexpr int WW = W_ + W_r;
  rfl.assign(5 * WW, 0.0);
  double angles[5];
  angles[0] = EMU[0];
  angles[1] = EMU[1];
  angles[2] = EMU[2];
  angles[3] = EMU[3];
  angles[4] = u0;
  for (int k = 0; k < NS2; ++k) {
    double wavel = state.WL[k];
    double osa_dir[5] = {};
    OSA::FJX_OSA(wavel, wind, chlr, angles, osa_dir);
    for (int j = 0; j < 5; ++j) {
      rfl[j + 5 * k] = albedo[j];
    }
  }
}

} // namespace Column
} // namespace CloudJ

#endif // CLOUDJ_COLUMN_BUILDER_HPP
