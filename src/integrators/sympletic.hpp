#pragma once

#include "../acceleration-accumulator.hpp"
#include "integrator.hpp"
#include "../omp_utils.hpp"

namespace integrators {

// sympletic euler integrator it's more accurate than basic euler
// diffrience is that it updates velocity first then position using NEW velocity
// that's why it's more accurate
template <int DIM, typename Attributes> 
class Sympletic : public Integrator<DIM, Attributes> {
private:
  AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator;

public:
  Sympletic(AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator)
      : accelerationAccumulator(accelerationAccumulator) {}

  void step(Bodies<DIM, Attributes> &bodies, double dt) override {

    accelerationAccumulator.compute(bodies);
    size_t n = bodies.localSize();

    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      bodies.local(i).vel() += accelerationAccumulator.accel(i) * dt;
      bodies.local(i).pos() += bodies.local(i).vel() * dt;
    }
  }
};

} // namespace integrators
