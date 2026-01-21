#pragma once

#include "func.hpp"
#include <cmath>
#include <numeric>

namespace forces {

struct charge {
  double charge;
};

template <int DIM> struct coulomb : public force<DIM, charge> {

public:
  // Vacuum permittivity
  static constexpr double k = 8.987551785972e9;

  Vec<DIM> operator()(const body<DIM, charge> &subjectBody,
                      const body<DIM, charge> &exertingBody) const override {

    constexpr double eps = 1e-12;

    Vec<DIM> displacement;
    std::transform(exertingBody.pos().begin(), exertingBody.pos().end(),
                   subjectBody.pos().begin(), displacement.begin(),
                   std::minus{});
    double dist2 = std::inner_product(displacement.begin(), displacement.end(),
                                      displacement.begin(), eps);
    double coeff = -k * subjectBody.attributes().charge *
                   exertingBody.attributes().charge /
                   (dist2 * std::sqrt(dist2));

    std::transform(displacement.begin(), displacement.end(),
                   displacement.begin(), [coeff](auto r) { return coeff * r; });

    return displacement;
  }
};

} // namespace forces
