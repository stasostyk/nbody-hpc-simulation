#pragma once

#include "../acceleration-accumulator.hpp"
#include "integrator.hpp"

namespace integrators {

// basic euler integrator it's the simplest but it can drift from accuracy
// easily and overshhot or undershhot aproximations it works by updating the
// position using OLD velocity (wgy it's not accurate) then updating velocity
// using acceleration
template <int DIM> 
class Euler : Integrator<DIM> {
private:
  const AccelerationAccumulator<DIM> accelerationAccumulator;

public:
  Euler(AccelerationAccumulator<DIM> accelerationAccumulator)
      : accelerationAccumulator(accelerationAccumulator) {}

  void step(std::vector<body<DIM>> &bodies, double dt) override {

    accelerationAccumulator.compute(bodies);
    size_t n = bodies.size();
    for (int i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies[i].position[d] += bodies[i].velocity[d] * dt;
      }
    }

    for (int i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies[i].velocity[d] += accelerationAccumulator.accel(i)[d] * dt;
      }
    }
  }
};

} // namespace integrators
