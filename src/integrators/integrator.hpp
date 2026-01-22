#pragma once

#include "../body.hpp"

namespace integrators {

template <int DIM, typename Attributes> 
class Integrator {
public:
  virtual void step(Bodies<DIM, Attributes> &bodies, double dt) = 0;
};

} // namespace integrators
