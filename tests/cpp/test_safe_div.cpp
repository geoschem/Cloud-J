/**
 * Property-based tests for CloudJ::SAFE_DIV
 *
 * Property 10: SAFE_DIV Never Produces NaN or Infinity
 *   For any double-precision numerator and denominator values and for any
 *   finite alt_nan, alt_overflow, and alt_underflow values, SAFE_DIV SHALL
 *   return a finite value (never NaN, never infinity).
 *
 * Property 11: SAFE_DIV Correctness for Safe Division
 *   For any numerator and denominator where the division does not overflow,
 *   underflow, or produce NaN (i.e., denom != 0, exponent difference within
 *   bounds), SAFE_DIV SHALL return the exact value numer/denom.
 *
 * Validates: Requirements 10.1, 10.2, 10.3, 10.4
 */

#include <cloudj/safe_div.hpp>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <limits>

// Simple xorshift64 PRNG for property-based testing
struct PRNG {
    uint64_t state;

    explicit PRNG(uint64_t seed) : state(seed ? seed : 1) {}

    uint64_t next() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }

    // Generate a random double in full double range (including subnormals, large values)
    double random_double() {
        uint64_t bits = next();
        double d;
        std::memcpy(&d, &bits, sizeof(d));
        // If NaN or Inf is generated, map to a finite value for input generation
        if (!std::isfinite(d)) {
            // Map to a large but finite value
            d = (bits & 1) ? DBL_MAX : -DBL_MAX;
        }
        return d;
    }

    // Generate a finite double (guaranteed)
    double random_finite_double() {
        double d;
        do {
            d = random_double();
        } while (!std::isfinite(d));
        return d;
    }

    // Generate double in a specific range [lo, hi]
    double random_in_range(double lo, double hi) {
        uint64_t bits = next();
        double t = static_cast<double>(bits) / static_cast<double>(UINT64_MAX);
        return lo + t * (hi - lo);
    }

    // Generate a "safe" pair: both non-zero, exponent difference within bounds
    void random_safe_pair(double& numer, double& denom) {
        // Generate values where division is safe:
        // Both non-zero, and exponent difference is well within max/min exponent bounds
        int max_exp = std::numeric_limits<double>::max_exponent;
        int min_exp = std::numeric_limits<double>::min_exponent;
        int safe_margin = 100; // keep well within bounds

        // Pick random exponents within a safe band
        int exp_n = static_cast<int>(next() % (max_exp - min_exp - 2 * safe_margin)) + min_exp + safe_margin;
        int exp_d = static_cast<int>(next() % (max_exp - min_exp - 2 * safe_margin)) + min_exp + safe_margin;

        // Ensure exponent difference is within bounds
        while (exp_n - exp_d >= max_exp || exp_n - exp_d <= min_exp) {
            exp_d = static_cast<int>(next() % (max_exp - min_exp - 2 * safe_margin)) + min_exp + safe_margin;
        }

        // Build doubles with these exponents: mantissa in [1,2)
        double mantissa_n = 1.0 + random_in_range(0.0, 1.0);
        double mantissa_d = 1.0 + random_in_range(0.0, 1.0);

        numer = std::ldexp(mantissa_n, exp_n);
        denom = std::ldexp(mantissa_d, exp_d);

        // Random sign
        if (next() & 1) numer = -numer;
        if (next() & 1) denom = -denom;
    }
};

static int failures = 0;

static void check(bool condition, const char* msg, int line) {
    if (!condition) {
        std::cerr << "FAIL (line " << line << "): " << msg << "\n";
        ++failures;
    }
}

#define CHECK(cond, msg) check((cond), (msg), __LINE__)

// ============================================================================
// Property 10: SAFE_DIV Never Produces NaN or Infinity
// Validates: Requirements 10.1, 10.2, 10.3, 10.4
// ============================================================================
static void test_property_10_never_nan_or_infinity() {
    std::cout << "Property 10: SAFE_DIV Never Produces NaN or Infinity\n";

    PRNG rng(42);

    // Fixed finite alternate values
    const double alt_nan = 0.0;
    const double alt_overflow = 999.0;
    const double alt_underflow = 0.0;

    // Test with random inputs
    const int N = 100000;
    for (int i = 0; i < N; ++i) {
        double numer = rng.random_finite_double();
        double denom = rng.random_finite_double();

        double result = CloudJ::SAFE_DIV(numer, denom, alt_nan, alt_overflow, alt_underflow);

        if (!std::isfinite(result)) {
            std::cerr << "  Counterexample: numer=" << numer
                      << " denom=" << denom
                      << " result=" << result << "\n";
            CHECK(false, "SAFE_DIV returned non-finite value for finite alternate values");
            return;
        }
    }

    // Targeted edge cases: 0/0, x/0, 0/x, huge/tiny, tiny/huge
    struct EdgeCase {
        double numer;
        double denom;
        const char* desc;
    };

    EdgeCase edge_cases[] = {
        {0.0, 0.0, "0/0"},
        {1.0, 0.0, "1/0"},
        {-1.0, 0.0, "-1/0"},
        {0.0, 1.0, "0/1"},
        {0.0, -1.0, "0/-1"},
        {DBL_MAX, DBL_MIN, "DBL_MAX/DBL_MIN"},
        {DBL_MIN, DBL_MAX, "DBL_MIN/DBL_MAX"},
        {-DBL_MAX, DBL_MIN, "-DBL_MAX/DBL_MIN"},
        {DBL_MIN, -DBL_MAX, "DBL_MIN/-DBL_MAX"},
        {DBL_MAX, 1.0, "DBL_MAX/1"},
        {1.0, DBL_MAX, "1/DBL_MAX"},
        {DBL_MAX, DBL_MAX, "DBL_MAX/DBL_MAX"},
        {DBL_MIN, DBL_MIN, "DBL_MIN/DBL_MIN"},
        {-DBL_MAX, -DBL_MAX, "-DBL_MAX/-DBL_MAX"},
        {1e308, 1e-308, "1e308/1e-308"},
        {1e-308, 1e308, "1e-308/1e308"},
        {-1e308, 1e-308, "-1e308/1e-308"},
        {5e-324, 5e-324, "smallest_subnormal/smallest_subnormal"},
        {DBL_MAX, -0.0, "DBL_MAX/-0"},
        {-0.0, -0.0, "-0/-0"},
        {-0.0, 0.0, "-0/0"},
        {0.0, -0.0, "0/-0"},
    };

    for (const auto& ec : edge_cases) {
        double result = CloudJ::SAFE_DIV(ec.numer, ec.denom, alt_nan, alt_overflow, alt_underflow);
        if (!std::isfinite(result)) {
            std::cerr << "  Counterexample [" << ec.desc << "]: result=" << result << "\n";
            CHECK(false, "SAFE_DIV returned non-finite value for edge case");
        }
    }

    // Test with varied finite alternate values
    for (int i = 0; i < 10000; ++i) {
        double numer = rng.random_finite_double();
        double denom = rng.random_finite_double();
        double a_nan = rng.random_finite_double();
        double a_ovf = rng.random_finite_double();
        double a_und = rng.random_finite_double();

        double result = CloudJ::SAFE_DIV(numer, denom, a_nan, a_ovf, a_und);

        if (!std::isfinite(result)) {
            std::cerr << "  Counterexample: numer=" << numer
                      << " denom=" << denom
                      << " alt_nan=" << a_nan
                      << " alt_overflow=" << a_ovf
                      << " alt_underflow=" << a_und
                      << " result=" << result << "\n";
            CHECK(false, "SAFE_DIV returned non-finite with random finite alternates");
            return;
        }
    }

    std::cout << "  PASSED (110000+ random + " 
              << (sizeof(edge_cases)/sizeof(edge_cases[0])) << " edge cases)\n";
}

// ============================================================================
// Property 11: SAFE_DIV Correctness for Safe Division
// Validates: Requirements 10.4
// ============================================================================
static void test_property_11_correctness_safe_division() {
    std::cout << "Property 11: SAFE_DIV Correctness for Safe Division\n";

    PRNG rng(123);

    const double alt_nan = -1.0;
    const double alt_overflow = -2.0;
    const double alt_underflow = -3.0;

    const int N = 100000;
    int tested = 0;

    for (int i = 0; i < N; ++i) {
        double numer, denom;
        rng.random_safe_pair(numer, denom);

        double result = CloudJ::SAFE_DIV(numer, denom, alt_nan, alt_overflow, alt_underflow);
        double expected = numer / denom;

        // For safe division, result should be exactly numer/denom
        if (result != expected) {
            std::cerr << "  Counterexample: numer=" << numer
                      << " denom=" << denom
                      << " expected=" << expected
                      << " got=" << result << "\n";
            CHECK(false, "SAFE_DIV did not return exact numer/denom for safe inputs");
            return;
        }
        ++tested;
    }

    // Also test specific known-safe cases
    struct SafeCase {
        double numer;
        double denom;
        const char* desc;
    };

    SafeCase safe_cases[] = {
        {1.0, 1.0, "1/1"},
        {-1.0, 1.0, "-1/1"},
        {1.0, -1.0, "1/-1"},
        {-1.0, -1.0, "-1/-1"},
        {0.0, 1.0, "0/1"},
        {0.0, -1.0, "0/-1"},
        {2.5, 0.5, "2.5/0.5"},
        {100.0, 3.0, "100/3"},
        {1e10, 1e5, "1e10/1e5"},
        {1e-10, 1e-5, "1e-10/1e-5"},
        {3.14159, 2.71828, "pi/e"},
        {-42.0, 7.0, "-42/7"},
    };

    for (const auto& sc : safe_cases) {
        double result = CloudJ::SAFE_DIV(sc.numer, sc.denom, alt_nan, alt_overflow, alt_underflow);
        double expected = sc.numer / sc.denom;
        if (result != expected) {
            std::cerr << "  Counterexample [" << sc.desc << "]: expected="
                      << expected << " got=" << result << "\n";
            CHECK(false, "SAFE_DIV did not return exact numer/denom for known-safe case");
        }
    }

    std::cout << "  PASSED (" << tested << " random safe pairs + "
              << (sizeof(safe_cases)/sizeof(safe_cases[0])) << " specific cases)\n";
}

// ============================================================================
// Additional boundary tests for Fortran behavior matching
// Validates: Requirements 10.1, 10.2, 10.3
// ============================================================================
static void test_boundary_behavior() {
    std::cout << "Boundary behavior tests (Fortran behavior matching)\n";

    // Req 10.1: 0/0 returns alt_nan
    {
        double result = CloudJ::SAFE_DIV(0.0, 0.0, 42.0, 99.0, -1.0);
        CHECK(result == 42.0, "0/0 should return alt_nan");
    }

    // Req 10.1: -0/0 is also 0/0
    {
        double result = CloudJ::SAFE_DIV(-0.0, 0.0, 42.0, 99.0, -1.0);
        CHECK(result == 42.0, "-0/0 should return alt_nan");
    }

    // Req 10.2: x/0 (non-zero numerator, zero denominator) returns alt_overflow
    {
        double result = CloudJ::SAFE_DIV(1.0, 0.0, 42.0, 99.0, -1.0);
        CHECK(result == 99.0, "1/0 should return alt_overflow");
    }
    {
        double result = CloudJ::SAFE_DIV(-5.0, 0.0, 42.0, 99.0, -1.0);
        CHECK(result == 99.0, "-5/0 should return alt_overflow");
    }
    {
        double result = CloudJ::SAFE_DIV(DBL_MAX, 0.0, 42.0, 99.0, -1.0);
        CHECK(result == 99.0, "DBL_MAX/0 should return alt_overflow");
    }

    // Req 10.2: overflow condition (exponent difference exceeds max_exponent)
    {
        double result = CloudJ::SAFE_DIV(DBL_MAX, DBL_MIN, 42.0, 99.0, -1.0);
        CHECK(result == 99.0, "DBL_MAX/DBL_MIN should return alt_overflow");
    }
    {
        double result = CloudJ::SAFE_DIV(1e308, 1e-308, 42.0, 99.0, -1.0);
        CHECK(result == 99.0, "1e308/1e-308 should return alt_overflow");
    }

    // Req 10.3: underflow condition (exponent difference below min_exponent)
    {
        double result = CloudJ::SAFE_DIV(DBL_MIN, DBL_MAX, 42.0, 99.0, -1.0);
        CHECK(result == -1.0, "DBL_MIN/DBL_MAX should return alt_underflow");
    }
    {
        double result = CloudJ::SAFE_DIV(1e-308, 1e308, 42.0, 99.0, -1.0);
        CHECK(result == -1.0, "1e-308/1e308 should return alt_underflow");
    }

    // Req 10.4: safe division returns exact result
    {
        double result = CloudJ::SAFE_DIV(6.0, 2.0, 42.0, 99.0, -1.0);
        CHECK(result == 3.0, "6/2 should return 3.0");
    }
    {
        double result = CloudJ::SAFE_DIV(1.0, 3.0, 42.0, 99.0, -1.0);
        CHECK(result == 1.0/3.0, "1/3 should return exact 1.0/3.0");
    }

    // Default argument behavior: alt_overflow defaults to NaN, alt_underflow defaults to 0.0
    // But since we are testing with finite alternates, those are covered above.

    std::cout << "  PASSED\n";
}

int main() {
    std::cout << "=== SAFE_DIV Property-Based Tests ===\n\n";

    test_property_10_never_nan_or_infinity();
    test_property_11_correctness_safe_division();
    test_boundary_behavior();

    std::cout << "\n=== Summary ===\n";
    if (failures == 0) {
        std::cout << "All tests PASSED\n";
        return 0;
    } else {
        std::cout << failures << " test(s) FAILED\n";
        return 1;
    }
}
