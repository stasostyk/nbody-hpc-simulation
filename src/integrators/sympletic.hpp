#pragma once

#include "../acceleration-accumulator.hpp"
#include "integrator.hpp"

namespace integrators {

// sympletic euler integrator it's more accurate than basic euler
// diffrience is that it updates velocity first then position using NEW velocity
// that's why it's more accurate
template <int DIM> 
class Sympletic : Integrator<DIM> {
private:
  const AccelerationAccumulator<DIM> accelerationAccumulator;

public:
  Sympletic(AccelerationAccumulator<DIM> accelerationAccumulator)
      : accelerationAccumulator(accelerationAccumulator) {}

  void step(std::vector<body<DIM>> &bodies, double dt) override {

    accelerationAccumulator.compute(bodies);
    size_t n = bodies.size();
    for (int i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies[i].velocity[d] += accelerationAccumulator.accel(i)[d] * dt;
        bodies[i].position[d] += bodies[i].velocity[d] * dt;
      }
    }
  }
};

} // namespace integrators
