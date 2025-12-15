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

  Vec<DIM> operator()(const body<DIM> &subjectBody,
                      const body<DIM> &exertingBody) const override {

    constexpr double eps = 1e-12;

    Vec<DIM> displacement;
    std::transform(exertingBody.pos().begin(), exertingBody.pos().end(),
                   subjectBody.pos().begin(), displacement.begin(), std::minus{});
    double dist2 = std::inner_product(displacement.begin(), displacement.end(),
                                      displacement.begin(), eps);
    double coeff =
        G * subjectBody.mass() * exertingBody.mass() / (dist2 * std::sqrt(dist2));

    std::transform(displacement.begin(), displacement.end(),
                   displacement.begin(), [coeff](auto r) { return coeff * r; });

    return displacement;
  }
};

} // namespace forces
