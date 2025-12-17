#pragma once

#include "../body.hpp"

namespace forces {

template <int DIM, typename Attributes> 
struct force {
public:
  virtual Vec<DIM> operator()(const body<DIM, Attributes> &subjectBody,
                              const body<DIM, Attributes> &exertingBody) const = 0;
};

} // namespace forces
