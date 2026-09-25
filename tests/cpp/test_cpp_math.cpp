#include <iostream>
#include <iomanip>
#include <cloudj/cross_sections.hpp>

int main() {
    double temp, t1, x1, t2, x2, t3, x3;
    int lqq;

    // Read inputs from stdin line-by-line matching Fortran test format
    while (std::cin >> temp >> t1 >> x1 >> t2 >> x2 >> t3 >> x3 >> lqq) {
        double res = CloudJ::CrossSections::interpolate(temp, t1, x1, t2, x2, t3, x3, lqq);
        std::cout << std::scientific << std::setprecision(15) << std::uppercase << res << "\n";
    }
    return 0;
}
