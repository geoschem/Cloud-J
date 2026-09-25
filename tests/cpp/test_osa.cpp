/// @file test_osa.cpp
/// @brief Property-based test for FJX_OSA: Ocean Surface Albedo bounds.
///
/// **Validates: Requirements 8.1, 8.4**
///
/// Property 13: Ocean Surface Albedo Bounds
/// For any wavelength in [200, 4000] nm, wind speed in [0, 30] m/s,
/// and chlorophyll-a in [0.01, 30] mg/m3, the computed ocean surface albedo
/// SHALL be in the range [0.0, 1.0].

#include <cloudj/osa.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

static int failures = 0;
static int total_tests = 0;

static void check_albedo_bounds(double wavel, double wind, double chlor,
                                const double* cangles, int num_angles) {
    double osa_dir[5] = {0.0};
    CloudJ::OSA::FJX_OSA(wavel, wind, chlor, cangles, osa_dir);

    for (int i = 0; i < 5; ++i) {
        ++total_tests;
        if (osa_dir[i] < 0.0 || osa_dir[i] > 1.0) {
            ++failures;
            std::cerr << "FAIL: OSA_dir[" << i << "] = " << std::setprecision(10)
                      << osa_dir[i] << " out of [0,1] range"
                      << " (wavel=" << wavel << " nm, wind=" << wind
                      << " m/s, chlor=" << chlor << " mg/m3"
                      << ", cos_angle=" << cangles[i] << ")\n";
        }
    }
}

int main() {
    // Test input ranges
    std::vector<double> wavelengths = {200, 300, 400, 500, 600, 700, 800, 1000, 2000, 3000, 4000};
    std::vector<double> wind_speeds = {0, 2, 5, 10, 15, 20, 25, 30};
    std::vector<double> chlorophylls = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0};

    // Cosine of zenith angles:
    // cos(0)=1.0, cos(30)=0.866, cos(60)=0.5, cos(80)=0.174, cos(89)=0.017
    std::vector<std::vector<double>> angle_sets = {
        {1.0, 0.866025, 0.5, 0.173648, 0.017452},       // standard set
        {0.017452, 0.017452, 0.017452, 0.017452, 0.017452}, // grazing incidence
        {0.5, 0.5, 0.5, 0.5, 0.5},                       // mid-latitude typical
    };

    std::cout << "=== Ocean Surface Albedo Bounds Test ===\n";
    std::cout << "Testing Property 13: all albedos in [0.0, 1.0]\n\n";

    // Exhaustive grid test across wavelengths, wind speeds, chlorophyll, angles
    for (double wavel : wavelengths) {
        for (double wind : wind_speeds) {
            for (double chlor : chlorophylls) {
                for (const auto& angles : angle_sets) {
                    check_albedo_bounds(wavel, wind, chlor, angles.data(), 5);
                }
            }
        }
    }

    // Edge case: very small angles (near grazing, cos ~ 0.001)
    {
        double grazing_angles[5] = {0.001, 0.005, 0.01, 0.05, 0.1};
        for (double wavel : wavelengths) {
            for (double wind : wind_speeds) {
                check_albedo_bounds(wavel, wind, 1.0, grazing_angles, 5);
            }
        }
    }

    // Edge case: cos(angle) = 0 (horizon)
    {
        double horizon_angles[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        for (double wavel : wavelengths) {
            check_albedo_bounds(wavel, 10.0, 1.0, horizon_angles, 5);
        }
    }

    // Edge case: typical mid-latitude ocean conditions
    // wavelength 500nm, wind 7 m/s, chlorophyll 0.3 mg/m3
    {
        double typical_angles[5] = {0.766, 0.643, 0.5, 0.342, 0.174};
        check_albedo_bounds(500.0, 7.0, 0.3, typical_angles, 5);
    }

    // Summary
    std::cout << "Total albedo values tested: " << total_tests << "\n";
    if (failures == 0) {
        std::cout << "PASSED: All albedo values are within [0.0, 1.0]\n";
    } else {
        std::cout << "FAILED: " << failures << " values out of range\n";
    }

    return (failures == 0) ? 0 : 1;
}
