#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <iostream>

// Vector in our DIM-dimensional space
template<int DIM>
struct Vec : public std::array<double, DIM> {
    using std::array<double, DIM>::array;  // Inherit constructors

    Vec(const std::array<double, DIM> &other) {
        for (int i = 0; i < DIM; i++) {
            (*this)[i] = other[i];
        }
    }

    Vec(double scalar) {
        for (int i = 0; i < DIM; i++) {
            (*this)[i] = scalar;
        }
    }

    Vec& operator=(double scalar) {
        std::fill(this->begin(), this->end(), scalar);
        return *this;
    }

    Vec& operator+=(const Vec& v) {
        std::transform(this->begin(), this->end(), v.begin(), this->begin(), std::plus{});
        return *this;
    }

    Vec& operator-=(const Vec& v) {
        std::transform(this->begin(), this->end(), v.begin(), this->begin(), std::minus{});
        return *this;
    }

    Vec& operator*=(const Vec& v) {
        std::transform(this->begin(), this->end(), v.begin(), this->begin(), std::multiplies{});
        return *this;
    }

    Vec operator+(const Vec& other) const {
        Vec result;
        for (size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] + other[i];
        }
        return result;
    }

    Vec operator-(const Vec& other) const {
        Vec result;
        for (size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] - other[i];
        }
        return result;
    }

    Vec& operator*=(double scalar) {
        for (size_t i = 0; i < DIM; ++i) {
            (*this)[i] *= scalar;
        }
        return *this;
    }

    Vec operator*(double scalar) const {
        Vec result = *this;
        result *= scalar;
        return result;
    }

    Vec& operator/=(double scalar) {
        for (size_t i = 0; i < DIM; ++i) {
            (*this)[i] /= scalar;
        }
        return *this;
    }

    Vec operator/(double scalar) const {
        Vec result = *this;
        result /= scalar;
        return result;
    }

    void print() const {
        std::cout << "(";
        for (int dimId = 0; dimId < DIM-1; dimId++)
            std::cout << (*this)[dimId] << ", ";
        std::cout << (*this)[DIM-1] << ")";
    }

    double normSquared() {
        double result = 0.;
        for (int i = 0; i < DIM; i++) {
            result += (*this)[i] * (*this)[i];
        }
        return result;
    }

    double norm() const {
        return sqrt(normSquared());
    }
};

// Non-member function for scalar * Vec
template<int DIM>
Vec<DIM> operator*(double scalar, const Vec<DIM>& v) {
    return v * scalar;  // Reuse the member function
}
