#pragma once

#include "../body.hpp"

namespace forces {

template <int DIM> 
struct force {
public:
  virtual Vec<DIM> operator()(const body<DIM> &subjectBody,
                              const body<DIM> &exertingBody) const = 0;
};

} // namespace forces
