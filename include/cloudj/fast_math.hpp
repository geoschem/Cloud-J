#ifndef CLOUDJ_FAST_MATH_HPP
#define CLOUDJ_FAST_MATH_HPP

#include <cmath>

// Thin wrapper around std::exp used throughout the radiative transfer code.
// Kept as a named function (rather than calling std::exp directly at each
// site) so that a platform-specific intrinsic (e.g. _mm256_exp_pd on AVX,
// or a GPU __expf path) can be swapped in at a single point without
// touching every call site.

namespace CloudJ {
namespace RadiativeSolver {

inline double exp_eval(double x) {
  return std::exp(x);
}

} // namespace RadiativeSolver
} // namespace CloudJ

#endif // CLOUDJ_FAST_MATH_HPP
