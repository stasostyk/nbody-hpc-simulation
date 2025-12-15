#pragma once

#include "../body.hpp"

namespace integrators {

template <int DIM> 
class Integrator {
public:
  virtual void step(bodies<DIM> &bodies, double dt) = 0;
};

} // namespace integrators
