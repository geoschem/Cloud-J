#ifndef CLOUDJ_SAFE_DIV_HPP
#define CLOUDJ_SAFE_DIV_HPP

#include <cmath>
#include <limits>

namespace CloudJ {

/// Perform "safe division" to prevent overflow, underflow, NaN, or infinity.
/// Returns an alternate value if the division cannot be performed safely.
///
/// Ported from SAFE_DIV in src/Core/cldj_error_mod.F90.
///
/// @param numer         Numerator
/// @param denom         Denominator
/// @param alt_nan       Value returned for 0/0 (NaN condition)
/// @param alt_overflow  Value returned on overflow (default: quiet NaN)
/// @param alt_underflow Value returned on underflow (default: 0.0)
inline double SAFE_DIV(double numer, double denom,
                       double alt_nan,
                       double alt_overflow = std::numeric_limits<double>::quiet_NaN(),
                       double alt_underflow = 0.0)
{
    // Case 1: 0/0 → return alt_nan
    if (numer == 0.0 && denom == 0.0) {
        return alt_nan;
    }

    // Extract binary exponents via std::frexp (equivalent to Fortran EXPONENT)
    int exp_n = 0;
    int exp_d = 0;
    std::frexp(numer, &exp_n);
    std::frexp(denom, &exp_d);

    // Case 2: Overflow — exponent difference exceeds max, or denom is zero
    if (exp_n - exp_d >= std::numeric_limits<double>::max_exponent
        || denom == 0.0) {
        return alt_overflow;
    }

    // Case 3: Underflow — exponent difference below min
    if (exp_n - exp_d <= std::numeric_limits<double>::min_exponent) {
        return alt_underflow;
    }

    // Case 4: Safe to divide
    return numer / denom;
}

} // namespace CloudJ

#endif // CLOUDJ_SAFE_DIV_HPP
