#ifndef CLOUDJ_RADIATIVE_SOLVER_HPP
#define CLOUDJ_RADIATIVE_SOLVER_HPP

#include <algorithm>
#include <array>
#include <cloudj/fast_math.hpp>
#include <cmath>
#include <experimental/mdspan.hpp>
#include <vector>

// Strictly opt-in OpenMP parallelization of the per-wavelength (18-bin) radiative
// solve inside a single MIESCT/OPMIE call. Disabled by default; only enabled
// when the CLOUDJ_USE_OPENMP compile definition is set (see CMakeLists.txt
// option CLOUDJ_USE_OPENMP). This is intended for callers that invoke
// cloud_jx() serially one column at a time and want to use idle CPU cores
// within a single call. Host models that already parallelize across columns
// should NOT enable this without benchmarking for oversubscription.
#if defined(CLOUDJ_USE_OPENMP)
#include <omp.h>
#endif

// Portable compiler loop-unrolling hint macros for standard compilers (GCC,
// Clang, Intel oneAPI)
#if defined(__clang__)
#define CLOUDJ_UNROLL_4 _Pragma("clang loop unroll_count(4)")
#elif defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER)
#define CLOUDJ_UNROLL_4 _Pragma("unroll(4)")
#elif defined(__GNUC__)
#define CLOUDJ_UNROLL_4 _Pragma("GCC unroll 4")
#else
#define CLOUDJ_UNROLL_4
#endif

namespace CloudJ {
namespace RadiativeSolver {

constexpr int M_ = 4;
constexpr int M2_ = 8;

namespace Photolysis {
constexpr int W_ = 18;
}

struct alignas(64) Workspace {
  // 2D buffers flat storage
  std::vector<double> a_data;  // size: M_ * nd
  std::vector<double> c_data;  // size: M_ * nd
  std::vector<double> h_data;  // size: M_ * nd
  std::vector<double> rr_data; // size: M_ * nd

  // 3D buffers flat storage
  std::vector<double> b_data;  // size: M_ * M_ * nd
  std::vector<double> aa_data; // size: M_ * M_ * nd
  std::vector<double> cc_data; // size: M_ * M_ * nd
  std::vector<double> dd_data; // size: M_ * M_ * nd

  // OPMIE integration buffers (hoisted out of OPMIE so the ~200KB of scratch
  // is allocated once per thread rather than malloc'd + zeroed on every call).
  // Only the shared OPMIE workspace carries these; the per-thread BLKSLV
  // scratch used inside the OpenMP MIESCT loop does not need them.
  std::vector<double> fj_data;     // size: nd * W_
  std::vector<double> fz_data;     // size: nd * W_
  std::vector<double> ztau_data;   // size: nd * W_
  std::vector<double> pomega_data; // size: M2_ * nd * W_

  // Resizes the OPMIE integration buffers to the current layer depth. Every
  // element of the [0,nd) x [active-k] region is written before it is read
  // (OPMIE setup fills even levels, the interpolation pass fills odd levels,
  // and BLKSLV fills fj), and inactive-k bins are never read, so no zeroing
  // is required — vector::resize only value-initialises newly-grown tail
  // elements and is a no-op in steady state.
  void resize_opmie(size_t nd) {
    size_t size_2d = nd * Photolysis::W_;
    fj_data.resize(size_2d);
    fz_data.resize(size_2d);
    ztau_data.resize(size_2d);
    pomega_data.resize(static_cast<size_t>(M2_) * size_2d);
  }

  // Resizes all vectors once to the required layer depth nd
  void resize(size_t nd) {
    size_t size_2d = M_ * nd;
    size_t size_3d = M_ * M_ * nd;

    // std::vector::resize is a no-op when the size is unchanged (steady
    // state), and only value-initialises newly-grown elements otherwise.
    // The active [0,nd) region of every buffer is fully written before it is
    // read: GEN_ID zeroes a/c/h/b/aa/cc in its init loop, and BLKSLV writes
    // dd/rr at every l. So there is no need to re-zero here, which is what
    // .assign() used to do on every call (and once per wavelength bin in the
    // OpenMP MIESCT path).
    a_data.resize(size_2d);
    c_data.resize(size_2d);
    h_data.resize(size_2d);
    rr_data.resize(size_2d);

    b_data.resize(size_3d);
    aa_data.resize(size_3d);
    cc_data.resize(size_3d);
    dd_data.resize(size_3d);
  }
};

// standard 8-stream Gauss points & weights from cldj_cmn_mod.F90
constexpr std::array<double, M_> EMU = {0.06943184420297, 0.33000947820757,
                                        0.66999052179243, 0.93056815579703};
constexpr std::array<double, M_> WT = {0.17392742256873, 0.32607257743127,
                                       0.32607257743127, 0.17392742256873};

using mdspan_2d =
    std::experimental::mdspan<const double,
                              std::experimental::dextents<size_t, 2>,
                              std::experimental::layout_left>;
using mdspan_1d =
    std::experimental::mdspan<const double,
                              std::experimental::dextents<size_t, 1>,
                              std::experimental::layout_left>;

using mdspan_3d_mut =
    std::experimental::mdspan<double, std::experimental::dextents<size_t, 3>,
                              std::experimental::layout_left>;
using mdspan_2d_mut =
    std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>,
                              std::experimental::layout_left>;
using mdspan_1d_mut =
    std::experimental::mdspan<double, std::experimental::dextents<size_t, 1>,
                              std::experimental::layout_left>;

// Static-extent variants for the fixed M_ = 4 / M2_ = 8 dimensions of the
// block-tridiagonal workspace: only the layer dimension `nd` stays dynamic.
// Index arithmetic (stride over l = M_*M_ etc.) folds to compile-time
// constants, enabling the same strength reduction gfortran gets from
// fixed-shape arrays. Offsets are arithmetically identical to the
// all-dynamic form, so results are bit-identical.
using mdspan_2d_mut_M = std::experimental::mdspan<
    double,
    std::experimental::extents<size_t, M_, std::experimental::dynamic_extent>,
    std::experimental::layout_left>;
using mdspan_2d_M2 = std::experimental::mdspan<
    const double,
    std::experimental::extents<size_t, M2_, std::experimental::dynamic_extent>,
    std::experimental::layout_left>;
using mdspan_3d_mut_MM = std::experimental::mdspan<
    double,
    std::experimental::extents<size_t, M_, M_,
                               std::experimental::dynamic_extent>,
    std::experimental::layout_left>;

/**
 * @brief Calculates ORDINARY Legendre functions of X
 * from P[0] = PL[0] = 1, P[1] = X, .... P[N-1] = PL[N-1]
 * Translates subroutine LEGND0 in cldj_fjx_sub_mod.F90.
 */
inline void LEGND0(double X, double PL[], int N) {
  PL[0] = 1.0;
  if (N > 1) {
    PL[1] = X;
    for (int i = 2; i < N; ++i) {
      double den = static_cast<double>(i);
      PL[i] = PL[i - 1] * X * (2.0 - 1.0 / den) - PL[i - 2] * (1.0 - 1.0 / den);
    }
  }
}

/**
 * @brief Generates coefficient matrices for the block tri-diagonal system.
 * Matches GEN_ID in cldj_fjx_sub_mod.F90.
 */
inline void GEN_ID(mdspan_2d_M2 pomega, // (M2_, N_)
                   mdspan_1d fz,        // (N_)
                   mdspan_1d ztau,      // (N_)
                   double zflux, const std::array<double, 5> &rfl,
                   const double pm[M_][M2_], const double pm0[M2_],
                   mdspan_3d_mut_MM b,  // (M_, M_, N_)
                   mdspan_3d_mut_MM aa, // (M_, M_, N_)
                   mdspan_3d_mut_MM cc, // (M_, M_, N_)
                   mdspan_2d_mut_M a,   // (M_, N_)
                   mdspan_2d_mut_M h,   // (M_, N_)
                   mdspan_2d_mut_M c,   // (M_, N_)
                   int nd) {
  // Local 4x4 matrix helpers
  double s[M_][M_] = {0};
  double t[M_][M_] = {0};
  double u[M_][M_] = {0};
  double v[M_][M_] = {0};
  double w[M_][M_] = {0};

  // Initialize outputs. The buffers are exactly M_*nd / M_*M_*nd elements, so
  // a contiguous std::fill covers precisely the [0,nd) region the strided
  // loop below used to walk — same values, vectorizable stores.
  std::fill(a.data_handle(), a.data_handle() + M_ * nd, 0.0);
  std::fill(c.data_handle(), c.data_handle() + M_ * nd, 0.0);
  std::fill(h.data_handle(), h.data_handle() + M_ * nd, 0.0);
  std::fill(b.data_handle(), b.data_handle() + M_ * M_ * nd, 0.0);
  std::fill(aa.data_handle(), aa.data_handle() + M_ * M_ * nd, 0.0);
  std::fill(cc.data_handle(), cc.data_handle() + M_ * M_ * nd, 0.0);

  // Upper boundary: 2nd-order terms
  int l1 = 0; // 0-based Fortran L=1
  int l2 = 1; // 0-based Fortran L=2

  // Hoist the loop-invariant pomega column loads: pomega is read-only but the
  // compiler cannot prove it does not alias the b/a/c/h outputs, so it would
  // otherwise re-issue these loads inside every i/j iteration.
  double po1[M2_], po2[M2_];
  for (int i = 0; i < M2_; ++i) {
    po1[i] = pomega(i, l1);
    po2[i] = pomega(i, l2);
  }

  for (int i = 0; i < M_; ++i) {
    double sum0 =
        po1[0] * pm[i][0] * pm0[0] + po1[2] * pm[i][2] * pm0[2] +
        po1[4] * pm[i][4] * pm0[4] + po1[6] * pm[i][6] * pm0[6];
    double sum2 =
        po2[0] * pm[i][0] * pm0[0] + po2[2] * pm[i][2] * pm0[2] +
        po2[4] * pm[i][4] * pm0[4] + po2[6] * pm[i][6] * pm0[6];
    double sum1 =
        po1[1] * pm[i][1] * pm0[1] + po1[3] * pm[i][3] * pm0[3] +
        po1[5] * pm[i][5] * pm0[5] + po1[7] * pm[i][7] * pm0[7];
    double sum3 =
        po2[1] * pm[i][1] * pm0[1] + po2[3] * pm[i][3] * pm0[3] +
        po2[5] * pm[i][5] * pm0[5] + po2[7] * pm[i][7] * pm0[7];
    h(i, l1) = 0.5 * (sum0 * fz(l1) + sum2 * fz(l2));
    a(i, l1) = 0.5 * (sum1 * fz(l1) + sum3 * fz(l2));
  }

  for (int i = 0; i < M_; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum0 = po1[0] * pm[i][0] * pm[j][0] +
                    po1[2] * pm[i][2] * pm[j][2] +
                    po1[4] * pm[i][4] * pm[j][4] +
                    po1[6] * pm[i][6] * pm[j][6];
      double sum2 = po2[0] * pm[i][0] * pm[j][0] +
                    po2[2] * pm[i][2] * pm[j][2] +
                    po2[4] * pm[i][4] * pm[j][4] +
                    po2[6] * pm[i][6] * pm[j][6];
      double sum1 = po1[1] * pm[i][1] * pm[j][1] +
                    po1[3] * pm[i][3] * pm[j][3] +
                    po1[5] * pm[i][5] * pm[j][5] +
                    po1[7] * pm[i][7] * pm[j][7];
      double sum3 = po2[1] * pm[i][1] * pm[j][1] +
                    po2[3] * pm[i][3] * pm[j][3] +
                    po2[5] * pm[i][5] * pm[j][5] +
                    po2[7] * pm[i][7] * pm[j][7];

      s[i][j] = -sum2 * WT[j];
      s[j][i] = -sum2 * WT[i];
      t[i][j] = -sum1 * WT[j];
      t[j][i] = -sum1 * WT[i];
      v[i][j] = -sum3 * WT[j];
      v[j][i] = -sum3 * WT[i];
      b(i, j, l1) = -0.5 * (sum0 + sum2) * WT[j];
      b(j, i, l1) = -0.5 * (sum0 + sum2) * WT[i];
    }
  }

  for (int i = 0; i < M_; ++i) {
    s[i][i] += 1.0;
    t[i][i] += 1.0;
    v[i][i] += 1.0;
    b(i, i, l1) += 1.0;

    c(i, l1) = s[i][0] * a(0, l1) / EMU[0] + s[i][1] * a(1, l1) / EMU[1] +
               s[i][2] * a(2, l1) / EMU[2] + s[i][3] * a(3, l1) / EMU[3];
  }

  for (int i = 0; i < M_; ++i) {
    for (int j = 0; j < M_; ++j) {
      w[j][i] = s[j][0] * t[0][i] / EMU[0] + s[j][1] * t[1][i] / EMU[1] +
                s[j][2] * t[2][i] / EMU[2] + s[j][3] * t[3][i] / EMU[3];
      u[j][i] = s[j][0] * v[0][i] / EMU[0] + s[j][1] * v[1][i] / EMU[1] +
                s[j][2] * v[2][i] / EMU[2] + s[j][3] * v[3][i] / EMU[3];
    }
  }

  double deltau = ztau(l2) - ztau(l1);
  double d2 = 0.25 * deltau;
  for (int i = 0; i < M_; ++i) {
    for (int j = 0; j < M_; ++j) {
      b(i, j, l1) += d2 * w[i][j];
      cc(i, j, l1) = d2 * u[i][j];
    }
    h(i, l1) += 2.0 * d2 * c(i, l1);
    a(i, l1) = 0.0;
  }

  for (int i = 0; i < M_; ++i) {
    double d1 = EMU[i] / deltau;
    b(i, i, l1) += d1;
    cc(i, i, l1) -= d1;
  }

  // Intermediate points: can be even or odd, A & C diagonal
  for (int ll = 1; ll <= nd - 2; ll += 2) {
    deltau = ztau(ll + 1) - ztau(ll - 1);
    const double p1 = pomega(1, ll), p3 = pomega(3, ll), p5 = pomega(5, ll),
                 p7 = pomega(7, ll);
    for (int i = 0; i < M_; ++i) {
      a(i, ll) = EMU[i] / deltau;
      c(i, ll) = -a(i, ll);
      h(i, ll) = fz(ll) * (p1 * pm[i][1] * pm0[1] + p3 * pm[i][3] * pm0[3] +
                           p5 * pm[i][5] * pm0[5] + p7 * pm[i][7] * pm0[7]);
    }
    for (int i = 0; i < M_; ++i) {
      for (int j = 0; j <= i; ++j) {
        double sum0 = p1 * pm[i][1] * pm[j][1] + p3 * pm[i][3] * pm[j][3] +
                      p5 * pm[i][5] * pm[j][5] + p7 * pm[i][7] * pm[j][7];
        b(i, j, ll) = -sum0 * WT[j];
        b(j, i, ll) = -sum0 * WT[i];
      }
    }
    for (int i = 0; i < M_; ++i) {
      b(i, i, ll) += 1.0;
    }
  }

  for (int ll = 2; ll <= nd - 3; ll += 2) {
    deltau = ztau(ll + 1) - ztau(ll - 1);
    const double p0 = pomega(0, ll), p2 = pomega(2, ll), p4 = pomega(4, ll),
                 p6 = pomega(6, ll);
    for (int i = 0; i < M_; ++i) {
      a(i, ll) = EMU[i] / deltau;
      c(i, ll) = -a(i, ll);
      h(i, ll) = fz(ll) * (p0 * pm[i][0] * pm0[0] + p2 * pm[i][2] * pm0[2] +
                           p4 * pm[i][4] * pm0[4] + p6 * pm[i][6] * pm0[6]);
    }
    for (int i = 0; i < M_; ++i) {
      for (int j = 0; j <= i; ++j) {
        double sum0 = p0 * pm[i][0] * pm[j][0] + p2 * pm[i][2] * pm[j][2] +
                      p4 * pm[i][4] * pm[j][4] + p6 * pm[i][6] * pm[j][6];
        b(i, j, ll) = -sum0 * WT[j];
        b(j, i, ll) = -sum0 * WT[i];
      }
    }
    for (int i = 0; i < M_; ++i) {
      b(i, i, ll) += 1.0;
    }
  }

  // Lower boundary: 2nd-order terms
  int l_last = nd - 1; // 0-based Fortran L=ND
  int l_prev = nd - 2; // 0-based Fortran L=ND-1

  double poL[M2_], poP[M2_];
  for (int i = 0; i < M2_; ++i) {
    poL[i] = pomega(i, l_last);
    poP[i] = pomega(i, l_prev);
  }

  for (int i = 0; i < M_; ++i) {
    double sum0 = poL[0] * pm[i][0] * pm0[0] + poL[2] * pm[i][2] * pm0[2] +
                  poL[4] * pm[i][4] * pm0[4] + poL[6] * pm[i][6] * pm0[6];
    double sum2 = poP[0] * pm[i][0] * pm0[0] + poP[2] * pm[i][2] * pm0[2] +
                  poP[4] * pm[i][4] * pm0[4] + poP[6] * pm[i][6] * pm0[6];
    double sum1 = poL[1] * pm[i][1] * pm0[1] + poL[3] * pm[i][3] * pm0[3] +
                  poL[5] * pm[i][5] * pm0[5] + poL[7] * pm[i][7] * pm0[7];
    double sum3 = poP[1] * pm[i][1] * pm0[1] + poP[3] * pm[i][3] * pm0[3] +
                  poP[5] * pm[i][5] * pm0[5] + poP[7] * pm[i][7] * pm0[7];
    h(i, l_last) = 0.5 * (sum0 * fz(l_last) + sum2 * fz(l_prev));
    a(i, l_last) = 0.5 * (sum1 * fz(l_last) + sum3 * fz(l_prev));
  }

  for (int i = 0; i < M_; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum0 = poL[0] * pm[i][0] * pm[j][0] + poL[2] * pm[i][2] * pm[j][2] +
                    poL[4] * pm[i][4] * pm[j][4] + poL[6] * pm[i][6] * pm[j][6];
      double sum2 = poP[0] * pm[i][0] * pm[j][0] + poP[2] * pm[i][2] * pm[j][2] +
                    poP[4] * pm[i][4] * pm[j][4] + poP[6] * pm[i][6] * pm[j][6];
      double sum1 = poL[1] * pm[i][1] * pm[j][1] + poL[3] * pm[i][3] * pm[j][3] +
                    poL[5] * pm[i][5] * pm[j][5] + poL[7] * pm[i][7] * pm[j][7];
      double sum3 = poP[1] * pm[i][1] * pm[j][1] + poP[3] * pm[i][3] * pm[j][3] +
                    poP[5] * pm[i][5] * pm[j][5] + poP[7] * pm[i][7] * pm[j][7];
      s[i][j] = -sum2 * WT[j];
      s[j][i] = -sum2 * WT[i];
      t[i][j] = -sum1 * WT[j];
      t[j][i] = -sum1 * WT[i];
      v[i][j] = -sum3 * WT[j];
      v[j][i] = -sum3 * WT[i];
      b(i, j, l_last) = -0.5 * (sum0 + sum2) * WT[j];
      b(j, i, l_last) = -0.5 * (sum0 + sum2) * WT[i];
    }
  }

  for (int i = 0; i < M_; ++i) {
    s[i][i] += 1.0;
    t[i][i] += 1.0;
    v[i][i] += 1.0;
    b(i, i, l_last) += 1.0;

    c(i, l_last) =
        s[i][0] * a(0, l_last) / EMU[0] + s[i][1] * a(1, l_last) / EMU[1] +
        s[i][2] * a(2, l_last) / EMU[2] + s[i][3] * a(3, l_last) / EMU[3];
  }

  for (int i = 0; i < M_; ++i) {
    for (int j = 0; j < M_; ++j) {
      w[j][i] = s[j][0] * t[0][i] / EMU[0] + s[j][1] * t[1][i] / EMU[1] +
                s[j][2] * t[2][i] / EMU[2] + s[j][3] * t[3][i] / EMU[3];
      u[j][i] = s[j][0] * v[0][i] / EMU[0] + s[j][1] * v[1][i] / EMU[1] +
                s[j][2] * v[2][i] / EMU[2] + s[j][3] * v[3][i] / EMU[3];
    }
  }

  deltau = ztau(l_last) - ztau(l_prev);
  d2 = 0.25 * deltau;
  double sumrfl = 0.0;
  for (int j = 0; j < M_; ++j) {
    sumrfl += rfl[j] * EMU[j] * WT[j];
  }
  double surfac = 4.0 / (1.0 + 2.0 * sumrfl);

  for (int i = 0; i < M_; ++i) {
    double d1 = EMU[i] / deltau;
    double sum0 = d1 + d2 * (w[i][0] + w[i][1] + w[i][2] + w[i][3]);
    for (int j = 0; j < M_; ++j) {
      aa(i, j, l_last) = -d2 * u[i][j];
      b(i, j, l_last) = b(i, j, l_last) + d2 * w[i][j] -
                        sum0 * surfac * rfl[j] * EMU[j] * WT[j];
    }
    h(i, l_last) = h(i, l_last) - 2.0 * d2 * c(i, l_last) +
                   sum0 * surfac * 0.25 * rfl[4] * zflux;
  }

  for (int i = 0; i < M_; ++i) {
    double d1 = EMU[i] / deltau;
    aa(i, i, l_last) += d1;
    b(i, i, l_last) += d1;
    c(i, l_last) = 0.0;
  }
}

// 4x4 LU solver helper matching manual BLKSLV algorithm (Highly Optimized
// Reciprocal Form)
inline void solve_lu_4x4(double E[M_][M_]) {
  double inv_E00 = 1.0 / E[0][0];
  E[1][0] *= inv_E00;
  E[1][1] = E[1][1] - E[1][0] * E[0][1];
  E[1][2] = E[1][2] - E[1][0] * E[0][2];
  E[1][3] = E[1][3] - E[1][0] * E[0][3];

  E[2][0] *= inv_E00;
  double inv_E11 = 1.0 / E[1][1];
  E[2][1] = (E[2][1] - E[2][0] * E[0][1]) * inv_E11;
  E[2][2] = E[2][2] - E[2][0] * E[0][2] - E[2][1] * E[1][2];
  E[2][3] = E[2][3] - E[2][0] * E[0][3] - E[2][1] * E[1][3];

  E[3][0] *= inv_E00;
  E[3][1] = (E[3][1] - E[3][0] * E[0][1]) * inv_E11;
  double inv_E22 = 1.0 / E[2][2];
  E[3][2] = (E[3][2] - E[3][0] * E[0][2] - E[3][1] * E[1][2]) * inv_E22;
  E[3][3] = E[3][3] - E[3][0] * E[0][3] - E[3][1] * E[1][3] - E[3][2] * E[2][3];

  // Invert L
  E[3][2] = -E[3][2];
  E[3][1] = -E[3][1] - E[3][2] * E[2][1];
  E[3][0] = -E[3][0] - E[3][1] * E[1][0] - E[3][2] * E[2][0];
  E[2][1] = -E[2][1];
  E[2][0] = -E[2][0] - E[2][1] * E[1][0];
  E[1][0] = -E[1][0];

  // Invert U (Using pre-calculated reciprocals for division-free speedups)
  E[3][3] = 1.0 / E[3][3];
  E[2][3] = -E[2][3] * E[3][3] * inv_E22;
  E[2][2] = inv_E22;
  E[1][3] = -(E[1][2] * E[2][3] + E[1][3] * E[3][3]) * inv_E11;
  E[1][2] = -E[1][2] * E[2][2] * inv_E11;
  E[1][1] = inv_E11;
  E[0][3] =
      -(E[0][1] * E[1][3] + E[0][2] * E[2][3] + E[0][3] * E[3][3]) * inv_E00;
  E[0][2] = -(E[0][1] * E[1][2] + E[0][2] * E[2][2]) * inv_E00;
  E[0][1] = -E[0][1] * E[1][1] * inv_E00;
  E[0][0] = inv_E00;

  // Multiply U-inverse * L-inverse, storing result in E
  double temp[M_][M_];
  for (int i = 0; i < M_; ++i) {
    for (int j = 0; j < M_; ++j) {
      temp[i][j] = E[i][j];
    }
  }

  E[0][0] = temp[0][0] + temp[0][1] * temp[1][0] + temp[0][2] * temp[2][0] +
            temp[0][3] * temp[3][0];
  E[0][1] = temp[0][1] + temp[0][2] * temp[2][1] + temp[0][3] * temp[3][1];
  E[0][2] = temp[0][2] + temp[0][3] * temp[3][2];
  E[1][0] = temp[1][1] * temp[1][0] + temp[1][2] * temp[2][0] +
            temp[1][3] * temp[3][0];
  E[1][1] = temp[1][1] + temp[1][2] * temp[2][1] + temp[1][3] * temp[3][1];
  E[1][2] = temp[1][2] + temp[1][3] * temp[3][2];
  E[2][0] = temp[2][2] * temp[2][0] + temp[2][3] * temp[3][0];
  E[2][1] = temp[2][2] * temp[2][1] + temp[2][3] * temp[3][1];
  E[2][2] = temp[2][2] + temp[2][3] * temp[3][2];
  E[3][0] = temp[3][3] * temp[3][0];
  E[3][1] = temp[3][3] * temp[3][1];
  E[3][2] = temp[3][3] * temp[3][2];
}

#if defined(CLOUDJ_USE_PCR)
// High-performance 4x4 matrix multiplication helper
inline void mat_mult_4x4(const double X[M_][M_], const double Y[M_][M_],
                         double Z[M_][M_]) {
  CLOUDJ_UNROLL_4
  for (int i = 0; i < M_; ++i) {
    CLOUDJ_UNROLL_4
    for (int j = 0; j < M_; ++j) {
      Z[i][j] = X[i][0] * Y[0][j] + X[i][1] * Y[1][j] + X[i][2] * Y[2][j] +
                X[i][3] * Y[3][j];
    }
  }
}

// High-performance 4x4 matrix by 4-vector multiplication helper
inline void mat_vec_mult_4(const double X[M_][M_], const double V[M_],
                           double R[M_]) {
  CLOUDJ_UNROLL_4
  for (int i = 0; i < M_; ++i) {
    R[i] = X[i][0] * V[0] + X[i][1] * V[1] + X[i][2] * V[2] + X[i][3] * V[3];
  }
}

// Inverts a 4x4 matrix in-place using our optimized division-free LU solver
inline void invert_matrix_4x4_helper(const double B[M_][M_],
                                     double invB[M_][M_]) {
  CLOUDJ_UNROLL_4
  for (int i = 0; i < M_; ++i) {
    CLOUDJ_UNROLL_4
    for (int j = 0; j < M_; ++j) {
      invB[i][j] = B[i][j];
    }
  }
  solve_lu_4x4(invB);
}

inline void solve_pcr(mdspan_2d_mut fj, mdspan_2d_M2 pomega, mdspan_1d fz,
                      mdspan_1d ztau, double fsbot,
                      const std::array<double, 5> &rfl,
                      const double pm[M_][M2_], const double pm0[M2_],
                      double &fjtop, double &fjbot,
                      std::array<double, 5> &fibot, int nd, int k_idx,
                      Workspace &ws) {
  // Create temporary block-tridiagonal views to assemble the global system
  // coefficients
  mdspan_2d_mut_M a(ws.a_data.data(), M_, nd);
  mdspan_2d_mut_M c(ws.c_data.data(), M_, nd);
  mdspan_2d_mut_M h(ws.h_data.data(), M_, nd);
  mdspan_2d_mut_M rr(ws.rr_data.data(), M_, nd);

  mdspan_3d_mut_MM b(ws.b_data.data(), M_, M_, nd);
  mdspan_3d_mut_MM aa(ws.aa_data.data(), M_, M_, nd);
  mdspan_3d_mut_MM cc(ws.cc_data.data(), M_, M_, nd);
  mdspan_3d_mut_MM dd(ws.dd_data.data(), M_, M_, nd);

  // Generate block tri-diagonal system coefficients (a, b, cc, c)
  GEN_ID(pomega, fz, ztau, fsbot, rfl, pm, pm0, b, aa, cc, a, h, c, nd);

  // Setup persistent Parallel Cyclic Reduction coefficient scratchpad vectors
  // To ensure zero allocation at runtime, we can utilize the dd_data and b_data
  // workspace pools or allocate small temporary buffers since this is
  // CPU-forced testing.
  std::vector<double> A_pcr(M_ * M_ * nd, 0.0);
  std::vector<double> B_pcr(M_ * M_ * nd, 0.0);
  std::vector<double> C_pcr(M_ * M_ * nd, 0.0);
  std::vector<double> D_pcr(M_ * nd, 0.0);

  mdspan_3d_mut_MM A_v(A_pcr.data(), M_, M_, nd);
  mdspan_3d_mut_MM B_v(B_pcr.data(), M_, M_, nd);
  mdspan_3d_mut_MM C_v(C_pcr.data(), M_, M_, nd);
  mdspan_2d_mut_M D_v(D_pcr.data(), M_, nd);

  // Initialize PCR matrices
  for (int l = 0; l < nd; ++l) {
    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
      D_v(i, l) = h(i, l);
      // a is diagonal of size 4 representing lower diagonal
      A_v(i, i, l) = a(i, l);
      CLOUDJ_UNROLL_4
      for (int j = 0; j < M_; ++j) {
        B_v(i, j, l) = b(i, j, l);
      }
    }
  }

  // Upper block C: level 0 uses full cc matrix, levels 1..nd-2 use diagonal c
  CLOUDJ_UNROLL_4
  for (int i = 0; i < M_; ++i) {
    CLOUDJ_UNROLL_4
    for (int j = 0; j < M_; ++j) {
      C_v(i, j, 0) = cc(i, j, 0);
    }
  }
  for (int l = 1; l < nd; ++l) {
    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
      C_v(i, i, l) = c(i, l);
    }
  }

  // Boundary conditions adjustments for Thomas equivalence
  // Thomas boundary elements are already generated in GEN_ID. We adjust lower
  // boundary directly.
  int l_last = nd - 1;
  for (int i = 0; i < M_; ++i) {
    for (int j = 0; j < M_; ++j) {
      A_v(i, j, l_last) = aa(i, j, l_last);
    }
  }

  // Execute O(log N) stride-reduction stages
  int stages = std::ceil(std::log2(nd));
  for (int step = 0; step < stages; ++step) {
    int stride = 1 << step;

    std::vector<double> A_next(M_ * M_ * nd, 0.0);
    std::vector<double> B_next(M_ * M_ * nd, 0.0);
    std::vector<double> C_next(M_ * M_ * nd, 0.0);
    std::vector<double> D_next(M_ * nd, 0.0);

    mdspan_3d_mut_MM An(A_next.data(), M_, M_, nd);
    mdspan_3d_mut_MM Bn(B_next.data(), M_, M_, nd);
    mdspan_3d_mut_MM Cn(C_next.data(), M_, M_, nd);
    mdspan_2d_mut_M Dn(D_next.data(), M_, nd);

    for (int l = 0; l < nd; ++l) {
      double alpha[M_][M_] = {0};
      double beta[M_][M_] = {0};

      // 1. Calculate alpha = - A_l * B_{l-stride}^-1
      if (l - stride >= 0) {
        double Bl_left[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Bl_left[i][j] = B_v(i, j, l - stride);
          }
        }
        double invBl_left[M_][M_];
        invert_matrix_4x4_helper(Bl_left, invBl_left);

        double Al[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Al[i][j] = A_v(i, j, l);
          }
        }
        double temp[M_][M_];
        mat_mult_4x4(Al, invBl_left, temp);
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            alpha[i][j] = -temp[i][j];
          }
        }
      }

      // 2. Calculate beta = - C_l * B_{l+stride}^-1
      if (l + stride < nd) {
        double Bl_right[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Bl_right[i][j] = B_v(i, j, l + stride);
          }
        }
        double invBl_right[M_][M_];
        invert_matrix_4x4_helper(Bl_right, invBl_right);

        double Cl[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Cl[i][j] = C_v(i, j, l);
          }
        }
        double temp[M_][M_];
        mat_mult_4x4(Cl, invBl_right, temp);
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            beta[i][j] = -temp[i][j];
          }
        }
      }

      // 3. Compute new coefficients (An, Bn, Cn, Dn)
      // Bn = B_l + alpha * C_{l-stride} + beta * A_{l+stride}
      double term1[M_][M_] = {0};
      double term2[M_][M_] = {0};
      if (l - stride >= 0) {
        double Cl_left[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Cl_left[i][j] = C_v(i, j, l - stride);
          }
        }
        mat_mult_4x4(alpha, Cl_left, term1);
      }
      if (l + stride < nd) {
        double Al_right[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Al_right[i][j] = A_v(i, j, l + stride);
          }
        }
        mat_mult_4x4(beta, Al_right, term2);
      }

      CLOUDJ_UNROLL_4
      for (int i = 0; i < M_; ++i) {
        CLOUDJ_UNROLL_4
        for (int j = 0; j < M_; ++j) {
          Bn(i, j, l) = B_v(i, j, l) + term1[i][j] + term2[i][j];
        }
      }

      // An = alpha * A_{l-stride}
      if (l - stride >= 0) {
        double Al_left[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Al_left[i][j] = A_v(i, j, l - stride);
          }
        }
        double temp[M_][M_];
        mat_mult_4x4(alpha, Al_left, temp);
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            An(i, j, l) = temp[i][j];
          }
        }
      }

      // Cn = beta * C_{l+stride}
      if (l + stride < nd) {
        double Cl_right[M_][M_];
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Cl_right[i][j] = C_v(i, j, l + stride);
          }
        }
        double temp[M_][M_];
        mat_mult_4x4(beta, Cl_right, temp);
        CLOUDJ_UNROLL_4
        for (int i = 0; i < M_; ++i) {
          CLOUDJ_UNROLL_4
          for (int j = 0; j < M_; ++j) {
            Cn(i, j, l) = temp[i][j];
          }
        }
      }

      // Dn = D_l + alpha * D_{l-stride} + beta * D_{l+stride}
      double termD1[M_] = {0};
      double termD2[M_] = {0};
      if (l - stride >= 0) {
        double Dl_left[M_];
        for (int i = 0; i < M_; ++i)
          Dl_left[i] = D_v(i, l - stride);
        mat_vec_mult_4(alpha, Dl_left, termD1);
      }
      if (l + stride < nd) {
        double Dl_right[M_];
        for (int i = 0; i < M_; ++i)
          Dl_right[i] = D_v(i, l + stride);
        mat_vec_mult_4(beta, Dl_right, termD2);
      }
      CLOUDJ_UNROLL_4
      for (int i = 0; i < M_; ++i) {
        Dn(i, l) = D_v(i, l) + termD1[i] + termD2[i];
      }
    }

    // Copy new coefficients back to PCR views for next stage
    for (int l = 0; l < nd; ++l) {
      CLOUDJ_UNROLL_4
      for (int i = 0; i < M_; ++i) {
        D_v(i, l) = Dn(i, l);
        CLOUDJ_UNROLL_4
        for (int j = 0; j < M_; ++j) {
          A_v(i, j, l) = An(i, j, l);
          B_v(i, j, l) = Bn(i, j, l);
          C_v(i, j, l) = Cn(i, j, l);
        }
      }
    }
  }

  // Final Stage: System is completely decoupled! Solve B'_l * X_l = D'_l
  for (int l = 0; l < nd; ++l) {
    double B_final[M_][M_];
    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
      CLOUDJ_UNROLL_4
      for (int j = 0; j < M_; ++j) {
        B_final[i][j] = B_v(i, j, l);
      }
    }
    double invB_final[M_][M_];
    invert_matrix_4x4_helper(B_final, invB_final);

    double D_final[M_];
    for (int i = 0; i < M_; ++i)
      D_final[i] = D_v(i, l);

    double X_final[M_];
    mat_vec_mult_4(invB_final, D_final, X_final);

    CLOUDJ_UNROLL_4
    for (int i = 0; i < M_; ++i) {
      rr(i, l) = X_final[i];
    }
  }

  // MEAN J & H (matching BLKSLV extraction exactly)
  for (int l = 0; l < nd; l += 2) {
    fj(l, k_idx) = rr(0, l) * WT[0] + rr(1, l) * WT[1] + rr(2, l) * WT[2] +
                   rr(3, l) * WT[3];
  }
  for (int l = 1; l < nd; l += 2) {
    fj(l, k_idx) = rr(0, l) * WT[0] * EMU[0] + rr(1, l) * WT[1] * EMU[1] +
                   rr(2, l) * WT[2] * EMU[2] + rr(3, l) * WT[3] * EMU[3];
  }

  // FJTOP diffuse flux out of top-of-atmosphere
  double sumt = rr(0, 0) * WT[0] * EMU[0] + rr(1, 0) * WT[1] * EMU[1] +
                rr(2, 0) * WT[2] * EMU[2] + rr(3, 0) * WT[3] * EMU[3];
  fjtop = 4.0 * sumt;

  // Surface diffuse flux integration
  double sumb = rr(0, l_last) * WT[0] * EMU[0] +
                rr(1, l_last) * WT[1] * EMU[1] +
                rr(2, l_last) * WT[2] * EMU[2] + rr(3, l_last) * WT[3] * EMU[3];

  double sumbr = rr(0, l_last) * WT[0] * EMU[0] * rfl[0] +
                 rr(1, l_last) * WT[1] * EMU[1] * rfl[1] +
                 rr(2, l_last) * WT[2] * EMU[2] * rfl[2] +
                 rr(3, l_last) * WT[3] * EMU[3] * rfl[3];

  double sumrf = WT[0] * EMU[0] * rfl[0] + WT[1] * EMU[1] * rfl[1] +
                 WT[2] * EMU[2] * rfl[2] + WT[3] * EMU[3] * rfl[3];

  double sumbx = (4.0 * sumbr + fsbot * rfl[4]) / (1.0 + 2.0 * sumrf);

  fjbot = 4.0 * sumb - sumbx;

  // FIBOT outputs: diffuse rays up/down
  fibot[4] = sumbx;
  for (int j = 0; j < 4; ++j) {
    fibot[j] = 2.0 * rr(j, l_last) - sumbx;
  }
}
#endif

/**
 * @brief Main radiative transfer tridiagonal solver logic.
 * Translates subroutine BLKSLV in cldj_fjx_sub_mod.F90.
 */
inline void BLKSLV(mdspan_2d_mut fj, // (N_, W_+W_r)
                   mdspan_2d_M2 pomega, // (M2_, N_) (from K-slice)
                   mdspan_1d fz,        // (N_) (from K-slice)
                   mdspan_1d ztau,      // (N_) (from K-slice)
                   double fsbot, const std::array<double, 5> &rfl,
                   const double pm[M_][M2_], const double pm0[M2_],
                   double &fjtop, double &fjbot, std::array<double, 5> &fibot,
                   int nd,
                   int k_idx,    // Current wavelength index
                   Workspace &ws // Persistent pre-allocated workspace reference
) {
#if defined(CLOUDJ_USE_PCR)
  // Redirect cleanly to our Parallel Cyclic Reduction solver backend
  solve_pcr(fj, pomega, fz, ztau, fsbot, rfl, pm, pm0, fjtop, fjbot, fibot, nd,
            k_idx, ws);
  return;
#endif

  // Create mdspan wrappers directly mapping over persistent workspace buffers
  // (zero allocation)
  mdspan_2d_mut_M a(ws.a_data.data(), M_, nd);
  mdspan_2d_mut_M c(ws.c_data.data(), M_, nd);
  mdspan_2d_mut_M h(ws.h_data.data(), M_, nd);
  mdspan_2d_mut_M rr(ws.rr_data.data(), M_, nd);

  mdspan_3d_mut_MM b(ws.b_data.data(), M_, M_, nd);
  mdspan_3d_mut_MM aa(ws.aa_data.data(), M_, M_, nd);
  mdspan_3d_mut_MM cc(ws.cc_data.data(), M_, M_, nd);
  mdspan_3d_mut_MM dd(ws.dd_data.data(), M_, M_, nd);

  // Generate block tri-diagonal system
  GEN_ID(pomega, fz, ztau, fsbot, rfl, pm, pm0, b, aa, cc, a, h, c, nd);




  // UPPER BOUNDARY L=1 (0-based: l=0)
  double E[M_][M_];
  for (int j = 0; j < M_; ++j) {
    for (int i = 0; i < M_; ++i) {
      E[i][j] = b(i, j, 0);
    }
  }

  solve_lu_4x4(E);

  for (int j = 0; j < M_; ++j) {
    for (int i = 0; i < M_; ++i) {
      dd(i, j, 0) = -E[i][0] * cc(0, j, 0) - E[i][1] * cc(1, j, 0) -
                    E[i][2] * cc(2, j, 0) - E[i][3] * cc(3, j, 0);
    }
    rr(j, 0) = E[j][0] * h(0, 0) + E[j][1] * h(1, 0) + E[j][2] * h(2, 0) +
               E[j][3] * h(3, 0);
  }

  // CONTINUE THROUGH ALL DEPTH POINTS l=1 TO nd-2 (Fortran L=2 TO ND-1)
  for (int l = 1; l < nd - 1; ++l) {
    for (int j = 0; j < M_; ++j) {
      for (int i = 0; i < M_; ++i) {
        b(i, j, l) += a(i, l) * dd(i, j, l - 1);
      }
      h(j, l) -= a(j, l) * rr(j, l - 1);
    }

    for (int j = 0; j < M_; ++j) {
      for (int i = 0; i < M_; ++i) {
        E[i][j] = b(i, j, l);
      }
    }


    solve_lu_4x4(E);

    for (int j = 0; j < M_; ++j) {
      for (int i = 0; i < M_; ++i) {
        dd(i, j, l) = -E[i][j] * c(j, l);
      }
      rr(j, l) = E[j][0] * h(0, l) + E[j][1] * h(1, l) + E[j][2] * h(2, l) +
                 E[j][3] * h(3, l);
    }
  }

  // FINAL DEPTH POINT: l=nd-1 (Fortran L=ND)
  int l_last = nd - 1;
  for (int j = 0; j < M_; ++j) {
    for (int i = 0; i < M_; ++i) {
      b(i, j, l_last) += aa(i, 0, l_last) * dd(0, j, l_last - 1) +
                         aa(i, 1, l_last) * dd(1, j, l_last - 1) +
                         aa(i, 2, l_last) * dd(2, j, l_last - 1) +
                         aa(i, 3, l_last) * dd(3, j, l_last - 1);
    }
    h(j, l_last) -= aa(j, 0, l_last) * rr(0, l_last - 1) +
                    aa(j, 1, l_last) * rr(1, l_last - 1) +
                    aa(j, 2, l_last) * rr(2, l_last - 1) +
                    aa(j, 3, l_last) * rr(3, l_last - 1);
  }

  for (int j = 0; j < M_; ++j) {
    for (int i = 0; i < M_; ++i) {
      E[i][j] = b(i, j, l_last);
    }
  }

  solve_lu_4x4(E);

  for (int j = 0; j < M_; ++j) {
    rr(j, l_last) = E[j][0] * h(0, l_last) + E[j][1] * h(1, l_last) +
                    E[j][2] * h(2, l_last) + E[j][3] * h(3, l_last);
  }

  // BACK SOLUTION
  for (int l = nd - 2; l >= 0; --l) {
    for (int j = 0; j < M_; ++j) {
      rr(j, l) += dd(j, 0, l) * rr(0, l + 1) + dd(j, 1, l) * rr(1, l + 1) +
                  dd(j, 2, l) * rr(2, l + 1) + dd(j, 3, l) * rr(3, l + 1);
    }
  }


  // MEAN J & H (Fortran 1-based level structure)
  // 0-based L indices: L=0, 2, 4... are odd levels (Fortran 1, 3, 5...)
  // L=1, 3, 5... are even levels (Fortran 2, 4, 6...)
  for (int l = 0; l < nd; l += 2) {
    fj(l, k_idx) = rr(0, l) * WT[0] + rr(1, l) * WT[1] + rr(2, l) * WT[2] +
                   rr(3, l) * WT[3];
  }
  for (int l = 1; l < nd; l += 2) {
    fj(l, k_idx) = rr(0, l) * WT[0] * EMU[0] + rr(1, l) * WT[1] * EMU[1] +
                   rr(2, l) * WT[2] * EMU[2] + rr(3, l) * WT[3] * EMU[3];
  }

  // FJTOP diffuse flux out of top-of-atmosphere
  double sumt = rr(0, 0) * WT[0] * EMU[0] + rr(1, 0) * WT[1] * EMU[1] +
                rr(2, 0) * WT[2] * EMU[2] + rr(3, 0) * WT[3] * EMU[3];
  fjtop = 4.0 * sumt;

  // Surface diffuse flux integration
  double sumb = rr(0, l_last) * WT[0] * EMU[0] +
                rr(1, l_last) * WT[1] * EMU[1] +
                rr(2, l_last) * WT[2] * EMU[2] + rr(3, l_last) * WT[3] * EMU[3];

  double sumbr = rr(0, l_last) * WT[0] * EMU[0] * rfl[0] +
                 rr(1, l_last) * WT[1] * EMU[1] * rfl[1] +
                 rr(2, l_last) * WT[2] * EMU[2] * rfl[2] +
                 rr(3, l_last) * WT[3] * EMU[3] * rfl[3];

  double sumrf = WT[0] * EMU[0] * rfl[0] + WT[1] * EMU[1] * rfl[1] +
                 WT[2] * EMU[2] * rfl[2] + WT[3] * EMU[3] * rfl[3];

  double sumbx = (4.0 * sumbr + fsbot * rfl[4]) / (1.0 + 2.0 * sumrf);

  fjbot = 4.0 * sumb - sumbx;

  // FIBOT outputs: diffuse rays up/down
  fibot[4] = sumbx;
  for (int j = 0; j < 4; ++j) {
    fibot[j] = 2.0 * rr(j, l_last) - sumbx;
  }

}

// Forward declaration of MIESCT
inline void MIESCT(mdspan_2d_mut fj,     // (N_, W_+W_r)
                   mdspan_1d_mut fjtop,  // (W_+W_r)
                   mdspan_1d_mut fjbot,  // (W_+W_r)
                   mdspan_2d_mut fibot,  // (5, W_+W_r)
                   mdspan_3d_mut pomega, // (M2_, N_, W_) (from Engine)
                   mdspan_2d_mut fz,     // (N_, W_) (from Engine)
                   mdspan_2d_mut ztau,   // (N_, W_) (from Engine)
                   mdspan_1d fsbot,      // (W_)
                   mdspan_2d rfl,        // (5, W_)
                   double u0, int nd,
                   const int *ldokr, // (W_+W_r) active-bin flags, Fortran LDOKR
                   Workspace &ws // Persistent pre-allocated workspace reference
);

/**
 * @brief Core light propagation and integration routine.
 * Translates subroutine OPMIE in cldj_fjx_sub_mod.F90.
 */
inline void OPMIE(mdspan_2d dtaux,       // (N_-1, W_)
                  mdspan_3d_mut pomegax, // (M2_, N_-1, W_)
                  double u0,
                  mdspan_2d rfl,                 // (5, W_)
                  mdspan_2d amf,                 // (N_, N_)
                  mdspan_1d amg,                 // (N_-1)
                  const std::vector<int> &jxtra, // (N_-1)
                  mdspan_2d_mut fjact,           // (N_-1, W_)
                  mdspan_1d_mut fjtop,           // (W_)
                  mdspan_1d_mut fjbot,           // (W_)
                  mdspan_2d_mut fibot,           // (5, W_)
                  mdspan_1d_mut fsbot,           // (W_)
                  mdspan_2d_mut fjflx,           // (N_-1, W_)
                  mdspan_2d_mut flxd,            // (N_-1, W_)
                  mdspan_1d_mut flxd0,           // (W_)
                  int lu,
                  const int *ldokr, // (W_+W_r) active-bin flags, Fortran LDOKR
                  double atau,      // geometric factor for sub-layer insertion
                  Workspace &ws // Persistent pre-allocated workspace reference
) {
  int l1u = lu + 1;
  int jaddto = 0;
  for (int l = 0; l < l1u; ++l) {
    jaddto += jxtra[l];
  }

  int nd = 2 * l1u + 2 * jaddto + 1;

  std::vector<int> l2lev(l1u + 1);
  l2lev[0] = 0; // 0-based
  for (int l = 1; l < l1u + 1; ++l) {
    int jx = (l - 1 < lu) ? jxtra[l - 1] : 0;
    l2lev[l] = l2lev[l - 1] + 1 + jx;
  }

  std::vector<double> ttau(l1u + 1, 0.0);
  std::vector<double> dtau1(l1u + 1, 0.0);
  std::vector<double> ftau(l1u + 1, 0.0);
  std::vector<double> pomega1_data(M2_ * (l1u + 1), 0.0);
  mdspan_2d_mut pomega1(pomega1_data.data(), M2_, l1u + 1);

  // MIESCT integration arrays: persistent workspace buffers (see
  // Workspace::resize_opmie) — no per-call allocation or zeroing.
  std::vector<double>& fj_data = ws.fj_data;
  std::vector<double>& fz_data = ws.fz_data;
  std::vector<double>& ztau_data = ws.ztau_data;
  std::vector<double>& pomega_data = ws.pomega_data;
  mdspan_2d_mut fj(fj_data.data(), nd, Photolysis::W_);
  mdspan_2d_mut fz(fz_data.data(), nd, Photolysis::W_);
  mdspan_2d_mut ztau(ztau_data.data(), nd, Photolysis::W_);
  mdspan_3d_mut pomega(pomega_data.data(), M2_, nd, Photolysis::W_);

  // Per-wavelength (k) setup: builds dtau1/ttau/ftau/pomega1 scratch and
  // scatters the results into the shared ztau/fz/pomega/flxd/flxd0 arrays.
  // ttau_l/dtau1_l/ftau_l/pomega1_l are fully overwritten at the start of
  // every call before being read (verified: dtau1 written for all l1u+1
  // entries; ttau written for all l1u+1 entries; ftau zeroed then
  // selectively written for all l1u+1 entries; pomega1 written for all
  // M2_*(l1u+1) entries) -- so they are safe per-k private scratch, but
  // MUST NOT be shared across threads (analogous to the Workspace race in
  // MIESCT). Everything else this lambda touches (dtaux, amg, amf, l2lev,
  // pomegax, and the shared ztau/fz/pomega/flxd/flxd0 outputs) is read or
  // written exclusively through the k-th column/slice, so different k
  // values never alias.
  auto opmie_setup_k = [&](int k, std::vector<double> &ttau_l,
                           std::vector<double> &dtau1_l,
                           std::vector<double> &ftau_l,
                           mdspan_2d_mut pomega1_l) {
    for (int l = 0; l < l1u; ++l) {
      dtau1_l[l] = dtaux(l, k) * amg(l);
    }
    dtau1_l[l1u] = 0.0;

    int ll0 = -1; // -1 represents no shadow boundary found
    for (int ll = 0; ll < l1u + 1; ++ll) {
      if (amf(ll, ll) <= 0.0) {
        ll0 = ll;
      }
    }

    for (int l = 0; l < l1u + 1; ++l)
      ftau_l[l] = 0.0;

    for (int ll = ll0 + 1; ll < l1u + 1; ++ll) {
      double xltau = 0.0;
      for (int ii = 0; ii < l1u; ++ii) {
        xltau += dtau1_l[ii] * amf(ii, ll);
      }
      if (xltau < 82.0) {
        ftau_l[ll] =
            exp_eval(-xltau); // Replaced std::exp with exp_eval dispatcher
      }
    }

    fsbot(k) = 0.0;
    if (ll0 == -1) {
      fsbot(k) = ftau_l[0] / amf(0, 0);
    }

    ttau_l[l1u] = 0.0;
    for (int l = l1u - 1; l >= 0; --l) {
      ttau_l[l] = ttau_l[l + 1] + dtaux(l, k) * (amg(l) * amg(l));
    }

    for (int i = 0; i < M2_; ++i) {
      pomega1_l(i, 0) = pomegax(i, 0, k);
      pomega1_l(i, l1u) = pomegax(i, l1u - 1, k);
    }
    for (int l = 1; l < l1u; ++l) {
      for (int i = 0; i < M2_; ++i) {
        pomega1_l(i, l) = (pomegax(i, l, k) * dtaux(l, k) +
                          pomegax(i, l - 1, k) * dtaux(l - 1, k)) /
                         (dtaux(l, k) + dtaux(l - 1, k));
      }
    }

    for (int l = 0; l < l1u + 1; ++l) {
      int l2 = l2lev[l];
      int lz = nd - 1 - 2 * l2;
      ztau(lz, k) = ttau_l[l];
      fz(lz, k) = ftau_l[l];
      for (int i = 0; i < M2_; ++i) {
        pomega(i, lz, k) = pomega1_l(i, l);
      }
    }

    // Geometric factor for sub-layer tau interpolation. The Fortran OPMIE
    // uses the module ATAU (the same value EXTRAL1 was called with, e.g.
    // 1.05); the source comment "ATAU = 1.05 / 0.005" documents the tuned
    // ATAU/ATAU0 *pair*, not a division. Using ATAU/ATAU0 here instead
    // overflows pow() for thick clouds and collapses the grid.

    for (int l = 0; l < l1u; ++l) {
      int l2 = l2lev[l];
      int lz = nd - 1 - 2 * l2;
      int l22 = l2lev[l + 1] - l2lev[l] - 1;

      if (l22 > 0) {
        double taubtm = ttau_l[l];
        double tautop = ttau_l[l + 1];
        double fbtm = ftau_l[l];
        double ftop = ftau_l[l + 1];
        double pombtm[M2_], pomtop[M2_];
        for (int i = 0; i < M2_; ++i) {
          pombtm[i] = pomega1_l(i, l);
          pomtop[i] = pomega1_l(i, l + 1);
        }

        double divt = 1.0 / (std::pow(atau, l22 + 1) - 1.0);

        for (int ll = 1; ll <= l22; ++ll) {
          int lzz = lz - 2 * ll;
          double sumt = (std::pow(atau, l22 + 1 - ll) - 1.0) * divt;
          ztau(lzz, k) = tautop + sumt * (taubtm - tautop);

          if (amf(0, 0) > 0.0) {
            double dtauext = (ztau(lzz, k) - tautop) / amg(l);
            fz(lzz, k) = ftop * exp_eval(-amf(l, l) *
                                         dtauext); // exp_eval dispatcher
          } else {
            double dtauext = (taubtm - ztau(lzz, k)) / amg(l);
            fz(lzz, k) = fbtm * exp_eval(-amf(l, l) *
                                         dtauext); // exp_eval dispatcher
          }

          for (int i = 0; i < M2_; ++i) {
            pomega(i, lzz, k) = pomtop[i] + sumt * (pombtm[i] - pomtop[i]);
          }
        }
      }
    }

    for (int l = ll0 + 1; l < l1u; ++l) {
      int l2 = l2lev[l];
      int lz = nd - 1 - 2 * l2;
      int l22 = l2lev[l + 1] - l2lev[l] - 1;

      if (l22 == 0) {
        double dtausca = dtau1_l[l] * pomegax(0, l, k);
        flxd(l, k) += 0.5 * (ftau_l[l] + ftau_l[l + 1]) * dtausca;

        double dtauabs = dtau1_l[l] - dtausca;
        double decay =
            exp_eval(-0.5 * dtauabs * amf(l, l)); // exp_eval dispatcher

        if (decay > 0.01) {
          flxd(l, k) += (ftau_l[l + 1] * (1.0 - decay) +
                         ftau_l[l] * (1.0 / decay - 1.0)) /
                        amf(l, l);
        } else {
          if (ll0 == -1) {
            flxd(l, k) += ftau_l[l + 1] / amf(l, l);
          } else {
            flxd(l, k) += ftau_l[l] / amf(l, l);
          }
        }
      } else {
        for (int ll = 0; ll <= l22; ++ll) {
          int lzz = lz - 2 * ll;
          double dtauext = (ztau(lzz, k) - ztau(lzz - 2, k)) / amg(l);
          double dtausca = dtauext * pomegax(0, l, k);
          flxd(l, k) += 0.5 * (fz(lzz, k) + fz(lzz - 2, k)) * dtausca;

          double dtauabs = dtauext - dtausca;
          double decay =
              exp_eval(-0.5 * dtauabs * amf(l, l)); // exp_eval dispatcher

          if (decay > 0.01) {
            flxd(l, k) += (fz(lzz - 2, k) * (1.0 - decay) +
                           fz(lzz, k) * (1.0 / decay - 1.0)) /
                          amf(l, l);
          } else {
            if (ll0 == -1) {
              flxd(l, k) += fz(lzz - 2, k) / amf(l, l);
            } else {
              flxd(l, k) += fz(lzz, k) / amf(l, l);
            }
          }
        }
      }
      flxd(l, k) *= amg(l);
    }

    for (int l = 0; l < l1u; ++l) {
      flxd0(k) += flxd(l, k);
    }

    for (int lz = 1; lz < nd - 1; lz += 2) {
      ztau(lz, k) = 0.5 * (ztau(lz - 1, k) + ztau(lz + 1, k));
      fz(lz, k) = std::sqrt(fz(lz - 1, k) * fz(lz + 1, k));
      for (int i = 0; i < M2_; ++i) {
        pomega(i, lz, k) = 0.5 * (pomega(i, lz - 1, k) + pomega(i, lz + 1, k));
      }
    }
  };

#if defined(CLOUDJ_USE_OPENMP)
  // Opt-in parallel path: each thread gets its own private ttau/dtau1/ftau/
  // pomega1 scratch (thread_local, resized on first touch per thread and
  // cheaply re-assigned every iteration -- no allocation once warmed up)
  // to avoid the data race described above. schedule(static) keeps
  // deterministic thread-to-iteration assignment; results are bit-identical
  // to the serial path regardless of schedule since k-slices are fully
  // independent.
#pragma omp parallel for schedule(static)
  for (int k = 0; k < Photolysis::W_; ++k) {
    if (ldokr[k] <= 0)
      continue; // Fortran: skip zero-solar bins (TROP-only NWBIN=8/12)
    thread_local std::vector<double> ttau_tls;
    thread_local std::vector<double> dtau1_tls;
    thread_local std::vector<double> ftau_tls;
    thread_local std::vector<double> pomega1_data_tls;
    ttau_tls.assign(l1u + 1, 0.0);
    dtau1_tls.assign(l1u + 1, 0.0);
    ftau_tls.assign(l1u + 1, 0.0);
    pomega1_data_tls.assign(static_cast<size_t>(M2_) * (l1u + 1), 0.0);
    mdspan_2d_mut pomega1_tls(pomega1_data_tls.data(), M2_, l1u + 1);

    opmie_setup_k(k, ttau_tls, dtau1_tls, ftau_tls, pomega1_tls);
  }
#else
  for (int k = 0; k < Photolysis::W_; ++k) {
    if (ldokr[k] <= 0)
      continue; // Fortran: skip zero-solar bins (TROP-only NWBIN=8/12)
    opmie_setup_k(k, ttau, dtau1, ftau, pomega1);
  }
#endif


  // Call MIESCT inside the workspace to orchestrate solving
  MIESCT(fj, fjtop, fjbot, fibot, pomega, fz, ztau, fsbot, rfl, u0, nd, ldokr,
         ws);

  // Post-processing (fjact/fjflx): reads/writes exclusively through the k-th
  // column of fj/fz/ztau/fjact/fjflx, and uses no shared mutable scratch, so
  // this loop is safe to parallelize directly.
  auto opmie_postproc_k = [&](int k) {
    for (int l = 0; l < l1u; ++l) {
      int lz0 = nd - 1 - 2 * l2lev[l + 1];
      int lz1 = nd - 1 - 2 * l2lev[l];

      double sumj =
          (4.0 * fj(lz0, k) + fz(lz0, k)) * (ztau(lz0 + 1, k) - ztau(lz0, k)) +
          (4.0 * fj(lz1, k) + fz(lz1, k)) * (ztau(lz1, k) - ztau(lz1 - 1, k));
      double sumt =
          ztau(lz0 + 1, k) - ztau(lz0, k) + ztau(lz1, k) - ztau(lz1 - 1, k);

      for (int lz = lz0 + 2; lz <= lz1 - 2; lz += 2) {
        sumj +=
            (4.0 * fj(lz, k) + fz(lz, k)) * (ztau(lz + 1, k) - ztau(lz - 1, k));
        sumt += ztau(lz + 1, k) - ztau(lz - 1, k);
      }
      fjact(l, k) = sumj / sumt;
    }

    for (int l = 1; l < l1u; ++l) {
      int lz = nd - 1 - 2 * l2lev[l];
      double fjflx0 =
          (ztau(lz + 1, k) - ztau(lz, k)) / (ztau(lz + 1, k) - ztau(lz - 1, k));
      fjflx(l - 1, k) =
          4.0 * (fj(lz - 1, k) * fjflx0 + fj(lz + 1, k) * (1.0 - fjflx0));
    }
  };

#if defined(CLOUDJ_USE_OPENMP)
#pragma omp parallel for schedule(static)
  for (int k = 0; k < Photolysis::W_; ++k) {
    if (ldokr[k] <= 0)
      continue; // Fortran: skipped bins keep their zeroed fjact/fjflx
    opmie_postproc_k(k);
  }
#else
  for (int k = 0; k < Photolysis::W_; ++k) {
    if (ldokr[k] <= 0)
      continue; // Fortran: skipped bins keep their zeroed fjact/fjflx
    opmie_postproc_k(k);
  }
#endif
}
inline void MIESCT(mdspan_2d_mut fj,     // (N_, W_+W_r)
                   mdspan_1d_mut fjtop,  // (W_+W_r)
                   mdspan_1d_mut fjbot,  // (W_+W_r)
                   mdspan_2d_mut fibot,  // (5, W_+W_r)
                   mdspan_3d_mut pomega, // (M2_, N_, W_) (from Engine)
                   mdspan_2d_mut fz,     // (N_, W_) (from Engine)
                   mdspan_2d_mut ztau,   // (N_, W_) (from Engine)
                   mdspan_1d fsbot,      // (W_)
                   mdspan_2d rfl,        // (5, W_)
                   double u0, int nd,
                   const int *ldokr, // (W_+W_r) active-bin flags, Fortran LDOKR
                   Workspace &ws // Persistent pre-allocated workspace reference
) {
  double pm[M_][M2_];
  double pm0[M2_];

  for (int i = 0; i < M_; ++i) {
    LEGND0(EMU[i], pm0, M2_);
    for (int im = 0; im < M2_; ++im) {
      pm[i][im] = pm0[im];
    }
  }

  // Note that U0 scattering does not change with altitude
  LEGND0(-u0, pm0, M2_);
  for (int im = 0; im < M2_; ++im) {
    pm0[im] = 0.25 * pm0[im];
  }

  // Body of the per-wavelength block solve, shared by both the serial and
  // OpenMP paths below. `local_ws` is the workspace to use for this k_idx:
  // the single shared `ws` in the serial case, or a private per-thread
  // workspace in the parallel case. Every array this touches (pomega_slice,
  // fz_slice, ztau_slice, fj, fjtop, fjbot, fibot) is indexed/offset by
  // k_idx, so no thread ever writes another thread's k_idx slice.
  auto miesct_solve_k = [&](int k_idx, Workspace &local_ws) {
    // Create views for the current wavelength slice to pass into BLKSLV
    mdspan_2d_M2 pomega_slice(pomega.data_handle() + k_idx * (M2_ * nd), M2_,
                              nd);
    mdspan_1d fz_slice(fz.data_handle() + k_idx * nd, nd);
    mdspan_1d ztau_slice(ztau.data_handle() + k_idx * nd, nd);

    // Map rfl array for the single wavelength
    std::array<double, 5> rfl_slice;
    for (int i = 0; i < 5; ++i) {
      rfl_slice[i] = rfl(i, k_idx);
    }

    std::array<double, 5> fibot_slice;
    double fjtop_val;
    double fjbot_val;

    BLKSLV(fj, pomega_slice, fz_slice, ztau_slice, fsbot(k_idx), rfl_slice, pm,
           pm0, fjtop_val, fjbot_val, fibot_slice, nd, k_idx, local_ws);

    fjtop(k_idx) = fjtop_val;
    fjbot(k_idx) = fjbot_val;
    for (int i = 0; i < 5; ++i) {
      fibot(i, k_idx) = fibot_slice[i];
    }
  };

#if defined(CLOUDJ_USE_OPENMP)
  // Opt-in parallel path: give every thread its OWN Workspace so the
  // internal scratch buffers (a_data, b_data, c_data, h_data, rr_data,
  // aa_data, cc_data, dd_data) are never shared/raced across threads. The
  // thread_local workspace is resized (cheap: .assign() on an
  // already-correctly-sized vector) on every iteration rather than only
  // once, since resize is idempotent and inexpensive, and this keeps the
  // logic simple and safe regardless of how OpenMP schedules iterations to
  // threads. schedule(static) gives deterministic, reproducible
  // thread-to-iteration assignment; results are bit-identical to the
  // serial path since each k_idx's solve is fully independent.
#pragma omp parallel for schedule(static)
  for (int k_idx = 0; k_idx < Photolysis::W_; ++k_idx) {
    if (ldokr[k_idx] <= 0)
      continue; // Fortran BLKSLV: skip zero-solar bins entirely
    thread_local Workspace tls_ws;
    tls_ws.resize(nd);
    miesct_solve_k(k_idx, tls_ws);
  }
#else
  // Default (safe) path: single shared workspace reused serially across all
  // 18 wavelength bins, exactly as before this change.
  for (int k_idx = 0; k_idx < Photolysis::W_; ++k_idx) {
    if (ldokr[k_idx] <= 0)
      continue; // Fortran BLKSLV: skip zero-solar bins entirely
    miesct_solve_k(k_idx, ws);
  }
#endif
}

} // namespace RadiativeSolver
} // namespace CloudJ

#endif // CLOUDJ_RADIATIVE_SOLVER_HPP
