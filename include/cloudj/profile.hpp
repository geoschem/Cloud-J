#ifndef CLOUDJ_PROFILE_HPP
#define CLOUDJ_PROFILE_HPP

#include <cloudj/column_builder.hpp>
#include <string>

namespace CloudJ {

// Input description for the high-level convenience API. It carries the same
// physical column information the standalone reference driver reads from
// atmos_PTClds.dat: surface conditions (month, latitude, surface pressure,
// broadband albedo, wind, chlorophyll), per-level climatology coefficients
// (eta-A/eta-B, temperature, relative humidity, aerosol loading), and the
// per-level cloud fields (fraction, liquid/ice water content). Everything the
// solver consumes (pressure/height edges, air/O3/CH4/H2O columns, cloud water
// paths and effective radii) is derived from these inputs by the shared
// column builder, so a profile fed to Engine::calculate_photolysis_rates
// produces exactly the column the reference driver would build from the same
// data.
class AtmosphericProfile {
public:
  AtmosphericProfile() = default;

  // Parse a reference-format column description (atmos_PTClds.dat layout).
  // Returns true on success.
  bool load_from_file(const std::string& path) {
    return Column::parse_atmos_ptclds(path, inputs_);
  }

  Column::Inputs& inputs() noexcept { return inputs_; }
  const Column::Inputs& inputs() const noexcept { return inputs_; }

private:
  Column::Inputs inputs_;
};

} // namespace CloudJ

#endif // CLOUDJ_PROFILE_HPP
