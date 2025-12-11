#pragma once

#include "func.hpp"
#include <cmath>
#include <numeric>

namespace forces {

template <int DIM> 
struct gravity : public force<DIM> {

public:
  // Gravitational constant
  static constexpr double G = 6.673e-11;

  Vec<DIM> operator()(const Vec<DIM> &subjectBodyPos, double subjectBodyMass,
                      const Vec<DIM> &exertingBodyPos,
                      double exertingBodyMass) const override {

    constexpr double eps = 1e-12;

    Vec<DIM> displacement;
    std::transform(exertingBodyPos.begin(), exertingBodyPos.end(),
                   subjectBodyPos.begin(), displacement.begin(), std::minus{});
    double dist2 = std::inner_product(displacement.begin(), displacement.end(),
                                      displacement.begin(), eps);
    double coeff =
        G * subjectBodyMass * exertingBodyMass / (dist2 * std::sqrt(dist2));

    std::transform(displacement.begin(), displacement.end(),
                   displacement.begin(), [coeff](auto r) { return coeff * r; });

    return displacement;
  }
};

} // namespace forces
