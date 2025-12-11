#pragma once

#include "../vec.hpp"

namespace forces {

template <int DIM> 
struct force {
public:
  virtual Vec<DIM> operator()(const Vec<DIM> &subjectBodyPos,
                              double subjectBodyMass,
                              const Vec<DIM> &exertingBodyPos,
                              double exertingBodyMass) const = 0;
};

} // namespace forces
