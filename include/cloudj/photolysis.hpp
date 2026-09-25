#ifndef CLOUDJ_PHOTOLYSIS_HPP
#define CLOUDJ_PHOTOLYSIS_HPP

#include <cloudj/cross_sections.hpp>
#include <cloudj/error.hpp>
#include <cmath>
#include <experimental/mdspan.hpp>
#include <sstream>
#include <string>
#include <vector>

namespace CloudJ {
namespace Photolysis {

constexpr int W_ = 18; // standard wavelengths
constexpr int S_ = W_; // v8.0: S_ == W_ == 18 (no Solar-J extension)

struct SpecData {
  int nw = W_;
  int ns = S_;
  int njx = 0;
  std::vector<std::string> titlejx;
  std::vector<char> sqq; // 'p' or 't' interpolation variable
  std::vector<int> lqq;  // number of points (1, 2, or 3)

  // interpolation temperatures/pressures per species, flat [njx * 3],
  // indexed as tqq[j * 3 + t]
  std::vector<double> tqq;

  // cross-sections, flat contiguous storage to avoid pointer-chasing in the
  // JRATET hot loop.
  std::vector<double> qo2; // [W_ * 3], indexed as qo2[k * 3 + t]
  std::vector<double> qo3; // [W_ * 3], indexed as qo3[k * 3 + t]
  std::vector<double> q1d; // [W_ * 3], indexed as q1d[k * 3 + t]
  std::vector<double>
      qqq; // [W_ * 3 * njx], indexed as qqq[(k * 3 + t) * njx + j]

  // Pre-computed reciprocal temperature span intervals
  std::vector<double>
      inv_t12; // [species] reciprocal of (tqq[species][1] - tqq[species][0])
  std::vector<double>
      inv_t23; // [species] reciprocal of (tqq[species][2] - tqq[species][1])

  inline double qo2_at(int k, int t) const { return qo2[k * 3 + t]; }
  inline double qo3_at(int k, int t) const { return qo3[k * 3 + t]; }
  inline double q1d_at(int k, int t) const { return q1d[k * 3 + t]; }
  inline double qqq_at(int k, int t, int j) const {
    return qqq[(k * 3 + t) * njx + j];
  }
  inline double tqq_at(int j, int t) const { return tqq[j * 3 + t]; }
};

/**
 * @brief Interpolates and accumulates photolysis rates (J-values) for a column.
 * Matches JRATET in cldj_fjx_sub_mod.F90.
 */
inline void
JRATET(const std::vector<double> &ppj, // pressure edges [lu + 1]
       const std::vector<double> &ttj, // mid-layer temperatures [lu + 1]
       std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>,
                                 std::experimental::layout_left>
           fff,                       // mean actinic fluxes [W_][lu]
       std::vector<double> &valjl,    // flat [lu * njxu], row-major [l][j]
       const SpecData &spec, int lu, int njxu, int &rc) {
  // Fortran JRATET (cldj_fjx_sub_mod.F90) errors out when the host model's
  // J-value array is too small: "CTM has not enough J-values dimensioned".
  if (njxu < spec.njx) {
    CLOUDJ_ERROR("CTM has not enough J-values dimensioned",
                 "JRATET in photolysis.hpp", rc);
    return;
  }
  // Flat, caller-reused buffer. Row-major layout valjl[l * njxu + j] matches
  // the previous nested valjl[l][j] exactly. We size + zero the whole buffer
  // (lu*njxu) so columns j in [spec.njx, njxu) stay 0, identical to the old
  // valjl.assign(lu, vector(njxu, 0.0)) behavior.
  valjl.assign(static_cast<size_t>(lu) * njxu, 0.0);

  // Per-species accumulator, reused across layers (allocated once per call
  // rather than once per layer; re-zeroed at the top of each iteration).
  std::vector<double> valj(spec.njx);

  for (int l = 0; l < lu; ++l) {
    double tt = ttj[l];
    double pp;
    if (l == 0) {
      pp = ppj[0];
    } else {
      pp = (ppj[l] + ppj[l + 1]) * 0.5;
    }

    // zero bin-11 below 100 hPa matching the O2 e-fold limit
    if (pp > 100.0) {
      fff(10, l) = 0.0; // 0-based index 10 corresponds to bin 11
    }

    std::fill(valj.begin(), valj.end(), 0.0);

    // Calculate O2, O3, and O3(1D) photolysis rates (reactions 0, 1, 2)
    for (int k = 0; k < W_; ++k) {
      double qo2tot = CrossSections::interpolate(
          tt, spec.tqq_at(0, 0), spec.qo2_at(k, 0), spec.tqq_at(0, 1),
          spec.qo2_at(k, 1), spec.tqq_at(0, 2), spec.qo2_at(k, 2), spec.lqq[0],
          spec.inv_t12[0], spec.inv_t23[0]);

      double qo3tot = CrossSections::interpolate(
          tt, spec.tqq_at(1, 0), spec.qo3_at(k, 0), spec.tqq_at(1, 1),
          spec.qo3_at(k, 1), spec.tqq_at(1, 2), spec.qo3_at(k, 2), spec.lqq[1],
          spec.inv_t12[1], spec.inv_t23[1]);

      double qo31dy = CrossSections::interpolate(
          tt, spec.tqq_at(2, 0), spec.q1d_at(k, 0), spec.tqq_at(2, 1),
          spec.q1d_at(k, 1), spec.tqq_at(2, 2), spec.q1d_at(k, 2), spec.lqq[2],
          spec.inv_t12[2], spec.inv_t23[2]);

      double qo31d = qo31dy * qo3tot;

      valj[0] += qo2tot * fff(k, l);
      valj[1] += qo3tot * fff(k, l);
      valj[2] += qo31d * fff(k, l);
    }

    // Calculate photolysis rates for reactions 4 to NJX (indices 3 to NJX-1)
    for (int j = 3; j < spec.njx; ++j) {
      double var = (spec.sqq[j] == 'p') ? pp : tt;
      for (int k = 0; k < W_; ++k) {
        double qqqt = CrossSections::interpolate(
            var, spec.tqq_at(j, 0), spec.qqq_at(k, 0, j), spec.tqq_at(j, 1),
            spec.qqq_at(k, 1, j), spec.tqq_at(j, 2), spec.qqq_at(k, 2, j),
            spec.lqq[j], spec.inv_t12[j], spec.inv_t23[j]);
        valj[j] += qqqt * fff(k, l);
      }
    }

    for (int j = 0; j < spec.njx; ++j) {
      valjl[l * njxu + j] = valj[j];
    }
  }
}

} // namespace Photolysis
} // namespace CloudJ

#endif // CLOUDJ_PHOTOLYSIS_HPP
