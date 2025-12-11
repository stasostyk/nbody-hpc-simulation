#pragma once

#include <algorithm>
#include <array>
#include <functional>

// Vector in our DIM-dimensional space
template<int DIM>
struct Vec : public std::array<double, DIM> {
    Vec& operator+=(const Vec& v) {
        std::transform(this->begin(), this->end(), v.begin(), this->begin(), std::plus{});
        return *this;
    }

    Vec& operator-=(const Vec& v) {
        std::transform(this->begin(), this->end(), v.begin(), this->begin(), std::minus{});
        return *this;
    }
};
