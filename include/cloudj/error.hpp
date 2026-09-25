#ifndef CLOUDJ_ERROR_HPP
#define CLOUDJ_ERROR_HPP

#include <string>
#include <stdexcept>
#include <iostream>

namespace CloudJ {

constexpr int CLDJ_SUCCESS = 0;
constexpr int CLDJ_FAILURE = -1;

class Error : public std::runtime_error {
public:
    explicit Error(const std::string& message) : std::runtime_error(message) {}
};

inline void CLOUDJ_ERROR(const std::string& message, const std::string& location, int& rc) {
    std::cerr << "CLOUDJ_ERROR: " << message << " --> LOCATION: " << location << "\n";
    rc = CLDJ_FAILURE;
}

inline void CLOUDJ_ERROR_STOP(const std::string& message, const std::string& location) {
    std::cerr << "CLOUDJ_ERROR_STOP: " << message << " --> LOCATION: " << location << "\n";
    throw Error(message + " at " + location);
}

} // namespace CloudJ

#endif // CLOUDJ_ERROR_HPP
