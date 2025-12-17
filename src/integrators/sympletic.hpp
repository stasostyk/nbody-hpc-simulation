#pragma once

#include "../acceleration-accumulator.hpp"
#include "integrator.hpp"

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

  void step(bodies<DIM, Attributes> &bodies, double dt) override {

    accelerationAccumulator.compute(bodies);
    size_t n = bodies.localSize();
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).vel()[d] += accelerationAccumulator.accel(i)[d] * dt;
        bodies.local(i).pos()[d] += bodies.local(i).vel()[d] * dt;
      }
    }
  }
};

} // namespace integrators
