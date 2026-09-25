// Cloud-J v8.0 C++ Standalone Driver
// Mirrors the Fortran cldj_standalone.F90 to produce identical output.

#ifndef MODEL_STANDALONE
#define MODEL_STANDALONE
#endif

#include <algorithm>
#include <chrono>
#include <cloudj/cloudj.hpp>
#include <cloudj/column_builder.hpp>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  std::string tables_dir = "tables";
  std::string output_path = "cpp_actual_output.txt";
  // bench_iters: number of in-memory repeated calls to engine.cloud_jx() to
  // time, excluding one-time init/atmosphere-setup cost. 0 = disabled
  // (default single-pass output-comparison mode).
  int bench_iters = 0;
  // init_tables_dir: optional directory of on-disk override tables passed to
  // engine.init() for the engine's internal spectral/aerosol/cloud tables
  // (FJX_spec.dat, FJX_scat-cld.dat, etc). Empty by default = embedded
  // (compiled-in) tables, matching current behavior. This is distinct from
  // `tables_dir` above, which locates the atmos_PTClds.dat profile file that
  // this driver reads directly and is unrelated to the engine's tables.
  std::string init_tables_dir = "";

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--output" && i + 1 < argc) {
      output_path = argv[++i];
    } else if (arg == "--bench-iters" && i + 1 < argc) {
      bench_iters = std::atoi(argv[++i]);
    } else if (arg == "--benchmark") {
      // Backward-compatible alias for benchmark/run_benchmark.py's
      // --benchmark-iters harness flag, which forwards a bare "--benchmark"
      // to this executable. Use a small default iteration count so the
      // harness still measures the corrected cloud_jx() throughput instead
      // of being a no-op.
      if (bench_iters == 0) bench_iters = 1000;
    } else if (arg == "--tables-dir" && i + 1 < argc) {
      // Only configures the engine's internal table loading (init_tables_dir).
      // We deliberately do NOT also repoint the existing `tables_dir` variable
      // (used for reading atmos_PTClds.dat) here: that file is not part of the
      // 10-table mapping this flag documents, and changing its resolution
      // based on this flag would be a surprising, undocumented side effect
      // for a flag named "tables-dir" in the context of engine table loading.
      init_tables_dir = argv[++i];
    }
  }

  // =====================================================================
  // Constants from state.hpp (MODEL_STANDALONE)
  // =====================================================================
  constexpr int L_  = CloudJ::CloudJState::L_;   // 57
  constexpr int L1_ = CloudJ::CloudJState::L1_;  // 58
  constexpr int L2_ = CloudJ::CloudJState::L2_;  // 59
  constexpr int LWEPAR = CloudJ::CloudJState::LWEPAR; // 34

  // =====================================================================
  // Engine initialization (mirrors Fortran INIT_CLDJ)
  // =====================================================================
  CloudJ::Engine engine;

  double ATAU_in    = 1.050;
  double ATAU0_in   = 0.005;
  double CLDCOR_in  = 0.33;
  int    NWBIN_in   = 18;
  int    LNRG_in    = 6;
  int    ATM0_in    = 1;
  int    CLDFLAG_in = 7;
  bool   use_h2o_uv_abs = true;

  int NJXX = 0;
  int rc = engine.init(
      true,       // amiroot
      L_,         // nlevels = 57
      LWEPAR,     // nlevels_with_cloud = 34
      CloudJ::JVN_,  // njxu
      ATAU_in, ATAU0_in, NWBIN_in,
      CLDFLAG_in, CLDCOR_in, LNRG_in, ATM0_in,
      use_h2o_uv_abs, NJXX,
      init_tables_dir);  // empty string by default = embedded tables

  if (rc != 0) {
    std::cerr << "Error: INIT_CLDJ failed with rc=" << rc << "\n";
    return 1;
  }

  // =====================================================================
  // Read atmosphere file: tables/atmos_PTClds.dat
  // =====================================================================
  // Parsing and column derivation live in the shared Column builder so the
  // library API and this reference driver stay in lock-step.
  CloudJ::Column::Inputs col_in;
  {
    std::string profile_filename = tables_dir + "/atmos_PTClds.dat";
    if (!CloudJ::Column::parse_atmos_ptclds(profile_filename, col_in)) {
      profile_filename = "./tables/atmos_PTClds.dat";
      if (!CloudJ::Column::parse_atmos_ptclds(profile_filename, col_in)) {
        std::cerr << "Error: Could not open profile file tables/atmos_PTClds.dat\n";
        return 1;
      }
    }
  }

  CloudJ::Column::Derived col;
  CloudJ::Column::build_column(col_in, engine.get_state(), col);

  const double* ALBEDO = col_in.albedo;
  const double WIND = col_in.wind;
  const double CHLR = col_in.chlr;

  // Aliases into the derived column, matching the names the solver calls below
  // use. CLF is mutated in place by CLOUD_JX and reset from CLF0 each call.
  double* const PPP = col.ppp.data();
  double* const ZZZ = col.zzz.data();
  double* const TTT = col.ttt.data();
  double* const DDD = col.ddd.data();
  double* const RRR = col.rrr.data();
  double* const OOO = col.ooo.data();
  double* const CCC = col.ccc.data();
  double* const HHH = col.hhh.data();
  double* const LWP = col.lwp.data();
  double* const IWP = col.iwp.data();
  double* const REFFL = col.reffl.data();
  double* const REFFI = col.reffi.data();
  double* const CLF = col.clf.data();
  const int* const CLDIW = col.cldiw.data();
  double* const AERSP = col.aersp.data();
  int* const NDXAER = col.ndxaer.data();
  const int LTOP = col.ltop;

  // Total spectral bins
  constexpr int WW = CloudJ::W_ + CloudJ::W_r;  // = 18 + 0 = 18

  // =====================================================================
  // In-memory benchmark mode: exercises the CORRECT, numerically-verified
  // engine.cloud_jx() path repeatedly, excluding the one-time init and
  // atmosphere-setup cost above (which a host model would only pay once
  // per column/timestep before calling cloud_jx). Reuses the same output
  // buffers across iterations, matching how a host model reuses buffers
  // across timesteps.
  // =====================================================================
  if (bench_iters > 0) {
    const CloudJ::CloudJState& state = engine.get_state();

    int NSZA = 30;  // mid-range, representative SZA
    double SZA = static_cast<double>(NSZA);
    int IRAN = 1;
    double SOLF = 1.0;
    double U0 = std::cos(SZA * CloudJ::CPI180);

    std::vector<double> RFL;
    CloudJ::Column::build_rfl(state, U0, ALBEDO, WIND, CHLR, RFL);

    // Disable diagnostic printing so I/O does not contaminate timing.
    bool LPRTJ = false;

    // Output buffers, reused (overwritten) across all iterations.
    double VALJXX[L_ * CloudJ::JVN_] = {};
    double SKPERD[(CloudJ::S_ + 2) * L1_] = {};
    double SWMSQ[6] = {};
    double OD18[L1_] = {};
    double WTQCA[CloudJ::NQD_] = {};
    int NICA = 0, JCOUNT = 0;
    bool LDARK = false;

    auto t_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < bench_iters; ++iter) {
      // CLOUD_JX mutates CLF in place, so reset it from the pristine copy each
      // call, matching how the normal SZA-scan loop resets it per iteration.
      CloudJ::Column::reset_cloud_fraction(col);
      rc = 0;
      engine.cloud_jx(
          U0, SZA, RFL.data(), SOLF, LPRTJ,
          PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
          LWP, IWP, REFFL, REFFI, CLF, CLDIW, CLDCOR_in,
          AERSP, NDXAER,
          L1_, CloudJ::AN_, CloudJ::JVN_,
          VALJXX, SKPERD, SWMSQ, OD18,
          IRAN, NICA, JCOUNT, LDARK,
          WTQCA, rc);
      if (rc != 0) {
        std::cerr << "Error: CLOUD_JX failed with rc=" << rc
                   << " at benchmark iteration " << iter << "\n";
        return 1;
      }
    }
    auto t_end = std::chrono::high_resolution_clock::now();

    double total_s =
        std::chrono::duration<double>(t_end - t_start).count();
    double per_call_ms = (total_s / static_cast<double>(bench_iters)) * 1000.0;
    double throughput = static_cast<double>(bench_iters) / total_s;

    std::cout << "BENCH: " << bench_iters << " iters, total="
              << std::fixed << std::setprecision(6) << total_s
              << " s, per-call=" << std::setprecision(6) << per_call_ms
              << " ms, throughput=" << std::setprecision(1) << throughput
              << " calls/s" << std::endl;

    return 0;
  }

  // =====================================================================
  // Print banner (matches Fortran)
  // =====================================================================
  std::cout << ">>>begin Cloud-J v8.0 Standalone" << std::endl;

  // =====================================================================
  // SZA scan loop (matches Fortran: 3 SZAs for test)
  // =====================================================================
  int SZAscan[3] = {0, 30, 60};

  for (int I = 0; I < 3; ++I) {
    int NSZA = SZAscan[I];
    double SZA = static_cast<double>(NSZA);

    // CLOUD_JX consumes CLF in place; reset from the pristine copy each iteration.
    CloudJ::Column::reset_cloud_fraction(col);

    int IRAN = 1;
    double SOLF = 1.0;
    double U0 = std::cos(SZA * CloudJ::CPI180);

    // Surface reflectivity: ocean OSA computed then overridden by the read-in
    // broadband albedo, exactly as the reference driver does.
    std::vector<double> RFL;
    const CloudJ::CloudJState& state = engine.get_state();
    CloudJ::Column::build_rfl(state, U0, ALBEDO, WIND, CHLR, RFL);

    // Set LPRTJ = true to trigger PHOTO_JX printing (matches Fortran)
    bool LPRTJ = true;

    // Output arrays
    double VALJXX[L_ * CloudJ::JVN_] = {};
    double SKPERD[(CloudJ::S_ + 2) * L1_] = {};
    double SWMSQ[6] = {};
    double OD18[L1_] = {};
    double WTQCA[CloudJ::NQD_] = {};
    int NICA = 0, JCOUNT = 0;
    bool LDARK = false;

    // Print SZA info (matches Fortran diagnostic)
    if (LPRTJ) {
      std::cout << std::fixed;
      std::cout << "SZA SOLF U0 albedo"
                << std::setw(8) << std::setprecision(3) << SZA
                << std::setw(8) << std::setprecision(5) << SOLF
                << std::setw(8) << std::setprecision(5) << U0
                << std::setw(8) << std::setprecision(5) << RFL[4 + 5 * (CloudJ::W_ - 1)]
                << std::endl;
      CloudJ::PhotoJX::JP_ATM0(PPP, TTT, DDD, OOO, ZZZ, L_);
      std::cout << " wvl  albedo u1:u4 & u0" << std::endl;
      for (int K = 0; K < CloudJ::NS2; ++K) {
        std::cout << std::setw(5) << (K + 1)
                  << std::fixed << std::setw(8) << std::setprecision(1) << state.WL[K];
        for (int J = 0; J < 5; ++J) {
          std::cout << std::fixed << std::setw(8) << std::setprecision(4) << RFL[J + 5 * K];
        }
        std::cout << std::endl;
      }
    }

    // Call CLOUD_JX
    rc = 0;
    engine.cloud_jx(
        U0, SZA, RFL.data(), SOLF, LPRTJ,
        PPP, ZZZ, TTT, HHH, DDD, RRR, OOO, CCC,
        LWP, IWP, REFFL, REFFI, CLF, CLDIW, CLDCOR_in,
        AERSP, NDXAER,
        L1_, CloudJ::AN_, CloudJ::JVN_,
        VALJXX, SKPERD, SWMSQ, OD18,
        IRAN, NICA, JCOUNT, LDARK,
        WTQCA, rc);

    if (rc != 0) {
      std::cerr << "Error: CLOUD_JX failed with rc=" << rc << "\n";
      return 1;
    }

    // The J-value table is printed from inside PHOTO_JX under LPRTJ, exactly
    // as the Fortran driver does: one table per SZA, from the first QCA.

  } // end SZA scan

  return 0;
}
