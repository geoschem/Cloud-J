// Unit tests for the high-level CloudJ::Engine library API.
//
// Checks are written with the CHECK() macro rather than assert() so they
// remain active in Release builds (NDEBUG); a dead assert silently turns a
// failing test into a passing one.

#include <cloudj/cloudj.hpp>
#include <cloudj/cross_sections.hpp>
#include <cloudj/error.hpp>
#include <cloudj/photolysis.hpp>
#include <cloudj/profile.hpp>
#include <cloudj/radiative_solver.hpp>
#include <cloudj/rates.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static int g_failures = 0;

static void check_expr(bool ok, const char *expr, const char *file, int line) {
  if (!ok) {
    std::cerr << "FAILED: " << expr << "  at " << file << ":" << line << "\n";
    ++g_failures;
  }
}
#define CHECK(cond) check_expr(static_cast<bool>(cond), #cond, __FILE__, __LINE__)

void test_profile_defaults() {
  // A default-constructed profile carries the standalone reference column
  // dimensions: L1_ per-level entries and LWEPAR cloud levels.
  CloudJ::AtmosphericProfile profile;
  CHECK(profile.inputs().etaa.size() ==
        static_cast<size_t>(CloudJ::CloudJState::L2_));
  CHECK(profile.inputs().cldfrw.size() ==
        static_cast<size_t>(CloudJ::CloudJState::LWEPAR));
}

void test_cross_section_interpolation() {
  // 1-point
  double x1 = CloudJ::CrossSections::interpolate(250.0, 200.0, 1.5, 0.0, 0.0,
                                                 0.0, 0.0, 1);
  CHECK(std::abs(x1 - 1.5) < 1e-10);

  // 2-point
  double x2 = CloudJ::CrossSections::interpolate(250.0, 200.0, 1.0, 300.0, 2.0,
                                                 0.0, 0.0, 2);
  CHECK(std::abs(x2 - 1.5) < 1e-10);

  // 2-point clamped lower
  double x2_low = CloudJ::CrossSections::interpolate(150.0, 200.0, 1.0, 300.0,
                                                     2.0, 0.0, 0.0, 2);
  CHECK(std::abs(x2_low - 1.0) < 1e-10);

  // 2-point clamped upper
  double x2_high = CloudJ::CrossSections::interpolate(350.0, 200.0, 1.0, 300.0,
                                                      2.0, 0.0, 0.0, 2);
  CHECK(std::abs(x2_high - 2.0) < 1e-10);

  // 3-point (lower half)
  double x3_low = CloudJ::CrossSections::interpolate(250.0, 200.0, 1.0, 300.0,
                                                     2.0, 400.0, 4.0, 3);
  CHECK(std::abs(x3_low - 1.5) < 1e-10);

  // 3-point (upper half)
  double x3_high = CloudJ::CrossSections::interpolate(350.0, 200.0, 1.0, 300.0,
                                                      2.0, 400.0, 4.0, 3);
  CHECK(std::abs(x3_high - 3.0) < 1e-10);
}

void test_radiative_solver_and_gauss_gauss() {
  // Check that we have EMU and WT parameters defined
  CHECK(CloudJ::RadiativeSolver::M_ == 4);
  CHECK(CloudJ::RadiativeSolver::M2_ == 8);
  CHECK(CloudJ::RadiativeSolver::EMU[0] > 0.0);
  CHECK(CloudJ::RadiativeSolver::WT[0] > 0.0);

  // Check manual LU block-diagonal tri-diagonal inversion math functions
  // compiles
  double E[CloudJ::RadiativeSolver::M_][CloudJ::RadiativeSolver::M_] = {
      {4.0, 1.0, 0.0, 0.0},
      {1.0, 4.0, 1.0, 0.0},
      {0.0, 1.0, 4.0, 1.0},
      {0.0, 0.0, 1.0, 4.0}};

  // Invert block using standard blockLU decomposition
  CloudJ::RadiativeSolver::solve_lu_4x4(E);

  // Verify E contains inverted matrix elements
  CHECK(E[0][0] > 0.0);
  CHECK(E[1][1] > 0.0);
}

// Parse one SZA block of the Fortran unit-7 diagnostic file (bin/fort.7),
// which the reference driver writes after CLOUD_JX. Each block has the
// header line "L Z p T [M] [O3] [H2O] [O1D] ..." followed by L_ data rows
// (top-down); column 8 is J1D = VALJXX(L,3), the O(1D) photolysis frequency
// of the FINAL cloud-weighted average. Returns the column indexed by 0-based
// layer (element 0 = level 1).
static std::vector<double> parse_fort7_j1d(const std::string &path, int block) {
  std::ifstream f(path);
  std::vector<double> out;
  if (!f.is_open()) return out;
  std::string line;
  int seen_blocks = -1;
  bool collecting = false;
  int rows = 0;
  while (std::getline(f, line)) {
    if (line.rfind("L Z p T", 0) == 0) {
      ++seen_blocks;
      collecting = (seen_blocks == block);
      rows = 0;
      continue;
    }
    if (!collecting) continue;
    std::istringstream ss(line);
    int lvl;
    if (!(ss >> lvl)) continue;
    std::vector<double> toks;
    double v;
    while (ss >> v) toks.push_back(v);
    if (toks.size() != 12) continue; // L + 12 data columns
    // toks[0..11] = ZKM, p, T, dAIR, dO3, dH2O, J1D, dO1D, POH, LCH4, HRuv, HRir
    if (lvl >= 1 && lvl <= CloudJ::CloudJState::L_) {
      if (static_cast<int>(out.size()) < CloudJ::CloudJState::L_)
        out.resize(CloudJ::CloudJState::L_, -1.0);
      out[lvl - 1] = toks[6];
      ++rows;
    }
    if (rows == CloudJ::CloudJState::L_) break;
  }
  return out;
}

void test_photolysis_and_orchestrator() {
  // The convenience API must run the real CLOUD_JX pipeline: for the same
  // physical column it must return exactly the J-values the reference
  // Fortran code produces, not a stand-in. The engine returns the FINAL
  // cloud-weighted ICA average, so the cross-language reference is the
  // Fortran unit-7 diagnostic (the driver's stdout J-table is the first-QCA
  // column printed from inside PHOTO_JX, not the average).
  CloudJ::Engine engine;
  int NJXX = 0;
  int rc = engine.init(true, CloudJ::CloudJState::L_,
                       CloudJ::CloudJState::LWEPAR, CloudJ::JVN_, 1.050,
                       0.005, 18, 7, 0.33, 6, 1, true, NJXX);
  CHECK(rc == CloudJ::CLDJ_SUCCESS);
  CHECK(NJXX > 0);

  CloudJ::AtmosphericProfile profile;
  std::string profile_path =
      std::string(CLOUDJ_TEST_TABLES_DIR) + "/atmos_PTClds.dat";
  CHECK(profile.load_from_file(profile_path));

  const int njx = engine.get_state().NJX;
  const int nlayers = CloudJ::CloudJState::L_;

  // Compute rates for SZA = 30 (the second block of the reference files).
  CloudJ::OutputRates rates = engine.calculate_photolysis_rates(profile, 30.0);
  CHECK(static_cast<int>(rates.j_values.size()) == nlayers);
  CHECK(static_cast<int>(rates.j_values[0].size()) == njx);

  // Rates must be physically sensible: finite, non-negative, and (daytime)
  // nonzero somewhere.
  bool any_positive = false;
  for (int l = 0; l < nlayers; ++l) {
    for (int j = 0; j < njx; ++j) {
      double v = rates.j_values[l][j];
      CHECK(std::isfinite(v));
      CHECK(v >= 0.0);
      if (v > 0.0) any_positive = true;
    }
  }
  CHECK(any_positive);

  // Compare the O(1D) column (species index 2) of the final average against
  // the Fortran unit-7 diagnostic for SZA=30. fort.7 uses E10.3 editing
  // (~4 sig figs), so a 1e-3 relative tolerance covers formatting only.
  std::vector<double> ref =
      parse_fort7_j1d(std::string(CLOUDJ_TEST_FORTRAN7), 1);
  CHECK(static_cast<int>(ref.size()) == nlayers);
  int mismatches = 0;
  for (int l = 0; l < nlayers; ++l) {
    double g = ref[l];
    if (g < 0.0) continue; // missing entry: skip, size check covers it
    double c = rates.j_values[l][2];
    double tol = 1e-3 * std::abs(g) + 1e-30;
    if (std::abs(g - c) > tol) {
      if (mismatches < 5) {
        std::cerr << "fort.7 mismatch L" << (l + 1) << ": ref=" << g
                  << " cpp=" << c << "\n";
      }
      ++mismatches;
    }
  }
  CHECK(mismatches == 0);

  // Dark conditions (SZA > 98) must give all-zero rates.
  CloudJ::OutputRates dark_rates =
      engine.calculate_photolysis_rates(profile, 100.0);
  CHECK(static_cast<int>(dark_rates.j_values.size()) == nlayers);
  for (int l = 0; l < nlayers; ++l)
    for (int j = 0; j < njx; ++j) CHECK(dark_rates.j_values[l][j] == 0.0);
}

void test_error_handling() {
  int rc = CloudJ::CLDJ_SUCCESS;
  CloudJ::CLOUDJ_ERROR("Warning test", "test_location", rc);
  CHECK(rc == CloudJ::CLDJ_FAILURE);

  bool threw = false;
  try {
    CloudJ::CLOUDJ_ERROR_STOP("Fatal error test", "test_location");
  } catch (const CloudJ::Error &e) {
    threw = true;
  }
  CHECK(threw);
}

void test_radiative_solver_workspace() {
  CloudJ::RadiativeSolver::Workspace ws;
  size_t nd = 10;
  ws.resize(nd);

  CHECK(ws.a_data.size() == CloudJ::RadiativeSolver::M_ * nd);
  CHECK(ws.c_data.size() == CloudJ::RadiativeSolver::M_ * nd);
  CHECK(ws.h_data.size() == CloudJ::RadiativeSolver::M_ * nd);
  CHECK(ws.rr_data.size() == CloudJ::RadiativeSolver::M_ * nd);

  CHECK(ws.b_data.size() ==
        CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);
  CHECK(ws.aa_data.size() ==
        CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);
  CHECK(ws.cc_data.size() ==
        CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);
  CHECK(ws.dd_data.size() ==
        CloudJ::RadiativeSolver::M_ * CloudJ::RadiativeSolver::M_ * nd);

  for (double val : ws.a_data) CHECK(val == 0.0);
  for (double val : ws.b_data) CHECK(val == 0.0);
}

int main() {
  std::cout << "Running standard library API unit tests...\n";
  test_profile_defaults();
  test_cross_section_interpolation();
  test_radiative_solver_and_gauss_gauss();
  test_photolysis_and_orchestrator();
  test_error_handling();
  test_radiative_solver_workspace();
  if (g_failures != 0) {
    std::cerr << g_failures << " check(s) FAILED\n";
    return 1;
  }
  std::cout << "All library API unit tests passed successfully!\n";
  return 0;
}
