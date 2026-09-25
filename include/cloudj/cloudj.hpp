#ifndef CLOUDJ_CLOUDJ_HPP
#define CLOUDJ_CLOUDJ_HPP

#include <array>
#include <cloudj/cloud_jx.hpp>
#include <cloudj/context.hpp>
#include <cloudj/cross_sections.hpp>
#include <cloudj/error.hpp>
#include <cloudj/init.hpp>
#include <cloudj/kokkos_backend.hpp>
#include <cloudj/osa.hpp>
#include <cloudj/photo_jx.hpp>
#include <cloudj/photolysis.hpp>
#include <cloudj/profile.hpp>
#include <cloudj/radiative_solver.hpp>
#include <cloudj/rates.hpp>
#include <cloudj/state.hpp>
#include <cmath>
#include <experimental/mdspan.hpp>
#include <string>
#include <vector>

namespace CloudJ {

using mdspan_2d_mut =
    std::experimental::mdspan<double, std::experimental::dextents<size_t, 2>,
                              std::experimental::layout_left>;

class Engine {
private:
  CloudJState state_;
  bool initialized_ = false;

public:
  Engine() = default;

  /**
   * @brief Computes photolysis rates (J-values) for a full atmospheric column.
   *
   * This is the high-level convenience entry point: it takes a physical column
   * description (AtmosphericProfile), derives every array the solver needs
   * using the same shared builder as the standalone reference driver, and runs
   * the real CLOUD_JX radiative-transfer pipeline. The returned J-values are
   * therefore identical to what the Fortran reference produces for the same
   * column, rather than a stand-in.
   *
   * init() must have been called first so the loaded tables and climatology
   * are available. The number of reactions per layer (j_values[l].size()) is
   * the engine's actual NJX, and the number of layers is the CTM layer count
   * (L_), matching the reference column.
   *
   * @param profile             Full column description (surface + per-level).
   * @param solar_zenith_angle  Solar zenith angle in degrees; > 98 deg is
   *                            treated as dark and returns zero rates.
   * @return                    J-values per layer per reaction.
   */
  OutputRates calculate_photolysis_rates(const AtmosphericProfile &profile,
                                         double solar_zenith_angle) {
    if (!initialized_) {
      throw Error(
          "calculate_photolysis_rates: engine not initialized; call init() "
          "first at Engine::calculate_photolysis_rates");
    }

    constexpr int L1 = CloudJState::L1_;
    constexpr int NLAYERS = CloudJState::L_;
    const int njx = state_.NJX;

    OutputRates rates;
    rates.j_values.assign(NLAYERS, std::vector<double>(njx, 0.0));

    // Build the derived column exactly as the standalone reference driver
    // would from the same inputs.
    Column::Derived col;
    Column::build_column(profile.inputs(), state_, col);

    // Dark shortcut mirrors CLOUD_JX/PHOTO_JX (SZA > 98 deg -> zero J-values),
    // so no solve is attempted when the column is in darkness.
    if (solar_zenith_angle > 98.0) {
      return rates;
    }

    double SZA = solar_zenith_angle;
    double U0 = std::cos(SZA * CPI180);
    double SOLF = 1.0;
    bool LPRTJ = false;
    int IRAN = 1;

    // Surface reflectivity for this geometry.
    std::vector<double> rfl;
    Column::build_rfl(state_, U0, profile.inputs().albedo,
                      profile.inputs().wind, profile.inputs().chlr, rfl);

    // Solver output buffers. VALJXX is laid out column-major [layer +
    // (L1U-1)*reaction]; the CLOUD_JX zeroing and accumulation assume the
    // (L1U-1) stride, so size it with the layer count and NJXU = JVN_.
    std::vector<double> valjxx(NLAYERS * JVN_, 0.0);
    std::vector<double> skperd((S_ + 2) * L1, 0.0);
    std::vector<double> swmsq(6, 0.0);
    std::vector<double> od18(L1, 0.0);
    std::vector<double> wtqca(NQD_, 0.0);
    int nica = 0, jcount = 0;
    bool ldark = false;
    int rc = CLDJ_SUCCESS;

    cloud_jx(U0, SZA, rfl.data(), SOLF, LPRTJ,
             col.ppp.data(), col.zzz.data(), col.ttt.data(), col.hhh.data(),
             col.ddd.data(), col.rrr.data(), col.ooo.data(), col.ccc.data(),
             col.lwp.data(), col.iwp.data(), col.reffl.data(),
             col.reffi.data(), col.clf.data(), col.cldiw.data(),
             state_.CLDCOR, col.aersp.data(), col.ndxaer.data(),
             L1, AN_, JVN_,
             valjxx.data(), skperd.data(), swmsq.data(), od18.data(),
             IRAN, nica, jcount, ldark, wtqca.data(), rc);

    if (rc != CLDJ_SUCCESS || ldark) {
      return rates; // j_values already zero-filled
    }

    // Unpack VALJXX (column-major: layer + NLAYERS*reaction) into the public
    // per-layer nested layout, using the engine's actual reaction count.
    for (int l = 0; l < NLAYERS; ++l) {
      for (int j = 0; j < njx; ++j) {
        rates.j_values[l][j] = valjxx[l + NLAYERS * j];
      }
    }
    return rates;
  }

  /**
   * @brief Initialize the Cloud-J engine by loading all lookup tables into
   * CloudJState. Wraps Init::INIT_CLDJ().
   *
   * @param amiroot    True if this is the root process (controls diagnostic output)
   * @param nlevels    Number of CTM levels
   * @param nlevels_with_cloud Number of CTM levels that include cloud data
   * @param njxu       Number of photolysis reactions requested by the caller
   * @param atau       Geometric factor for inserted sub-layers (default 1.120)
   * @param atau0      Minimum OD threshold for inserted layers (default 0.005)
   * @param nwbin      Number of wavelength bins to use (8, 12, or 18)
   * @param cldflag    Cloud overlap scheme flag (1-8, not 4)
   * @param cldcor     Cloud decorrelation parameter
   * @param lnrg       Max-overlap group selection (0, 3, or 6)
   * @param atm0       Atmosphere type selector for reference profiles
   * @param use_h2o_uv_abs Enable H2O UV absorption
   * @param njxx       Output: actual number of photolysis reactions available
   * @param datadir    Optional path to a directory containing on-disk override
   *                   tables (FJX_spec.dat, FJX_scat-cld.dat, etc. -- see
   *                   init.hpp for the full filename list). If empty (default)
   *                   or a specific file is not found there, falls back to the
   *                   embedded compiled-in tables. A message is printed to
   *                   stderr only when a disk file is actually used.
   * @return           CLDJ_SUCCESS (0) on success, non-zero on error
   */
  int init(bool amiroot, int nlevels, int nlevels_with_cloud,
           int njxu, double atau, double atau0, int nwbin,
           int cldflag, double cldcor, int lnrg, int atm0,
           bool use_h2o_uv_abs, int& njxx,
           const std::string& datadir = "") {
    int rc = CLDJ_SUCCESS;
    std::vector<std::string> titlejxx;

    Init::INIT_CLDJ(
        amiroot, datadir,
        nlevels, nlevels_with_cloud,
        titlejxx, njxu,
        atau, atau0, nwbin,
        cldflag, cldcor, lnrg, atm0,
        use_h2o_uv_abs,
        njxx, state_, rc);

    if (rc == CLDJ_SUCCESS) {
      initialized_ = true;
    }
    return rc;
  }

  /**
   * @brief Execute the full CLOUD_JX pipeline for a single atmospheric column.
   * Wraps the global CLOUD_JX() function using the Engine's internal CloudJState.
   *
   * @param u0         cos(solar zenith angle)
   * @param sza        Solar zenith angle (degrees)
   * @param rfl_flat   Surface reflectivity array [5 * (W_ + W_r)]
   * @param solf       Solar flux factor
   * @param lprtj      Enable diagnostic printing
   * @param ppp        Pressure edges [L1U+1]
   * @param zzz        Height edges [L1U+1] (cm)
   * @param ttt        Temperature per layer [L1U]
   * @param hhh        H2O column per layer [L1U]
   * @param ddd        Air density per layer [L1U]
   * @param rrr        Relative humidity per layer [L1U]
   * @param ooo        O3 column per layer [L1U]
   * @param ccc        CH4 column per layer [L1U]
   * @param lwp        Liquid water path per layer [L1U]
   * @param iwp        Ice water path per layer [L1U]
   * @param reffl      Liquid cloud effective radius [L1U]
   * @param reffi      Ice cloud effective radius [L1U]
   * @param cldf       Cloud fraction per layer [L1U]
   * @param cldiw      Cloud type index per layer [L1U]
   * @param cldcor_in  Cloud decorrelation parameter
   * @param aersp      Aerosol species paths [AN_ * L1U]
   * @param ndxaer     Aerosol type indices [AN_ * L1U]
   * @param l1u        Number of levels (L1U = L_ + 1)
   * @param anu        Number of aerosol types in layer (AN_)
   * @param njxu       Number of photolysis reactions
   * @param valjxx     Output: J-values [(L1U-1) * NJXU]
   * @param skperd     Output: heating rates [(S_+2) * L1U]
   * @param swmsq      Output: 6 solar flux diagnostics
   * @param od18       Output: per-layer OD at 18 bins [L1U]
   * @param iran       Random number seed index
   * @param nica       Output: number of ICAs computed
   * @param jcount     Output: number of PHOTO_JX calls made
   * @param ldark      Output: true if column is dark (SZA > limit)
   * @param wtqca      Output: quadrature weights [NQD_]
   * @param rc         Output: return code (0=success)
   * @param dir_sfc_flux   Optional output: direct surface flux [W_]
   * @param diff_sfc_flux  Optional output: diffuse surface flux [W_]
   * @param dep_flux       Optional output: deposition flux [W_]
   * @param diff_top_flux  Optional output: diffuse top-of-atmosphere flux [W_]
   */
  void cloud_jx(
      double u0, double sza,
      const double* rfl_flat, double solf,
      bool lprtj,
      const double* ppp, const double* zzz,
      const double* ttt, const double* hhh,
      const double* ddd, const double* rrr,
      const double* ooo, const double* ccc,
      const double* lwp, const double* iwp,
      const double* reffl, const double* reffi,
      const double* cldf, const int* cldiw,
      double cldcor_in,
      const double* aersp, const int* ndxaer,
      int l1u, int anu, int njxu,
      double* valjxx, double* skperd, double* swmsq, double* od18,
      int iran, int& nica, int& jcount, bool& ldark,
      double* wtqca,
      int& rc,
      double* dir_sfc_flux = nullptr,
      double* diff_sfc_flux = nullptr,
      double* dep_flux = nullptr,
      double* diff_top_flux = nullptr) {

    CLOUD_JX(u0, sza, rfl_flat, solf, lprtj,
             ppp, zzz, ttt, hhh, ddd, rrr, ooo, ccc,
             lwp, iwp, reffl, reffi, cldf, cldiw,
             cldcor_in,
             aersp, ndxaer,
             l1u, anu, njxu,
             valjxx, skperd, swmsq, od18,
             iran, nica, jcount, ldark,
             wtqca,
             state_, rc,
             dir_sfc_flux, diff_sfc_flux, dep_flux, diff_top_flux);
  }

  /**
   * @brief Returns a const reference to the internal CloudJState.
   */
  const CloudJState& get_state() const noexcept { return state_; }
};

} // namespace CloudJ

#endif // CLOUDJ_CLOUDJ_HPP
