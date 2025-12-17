#pragma once

#include "../acceleration-accumulator.hpp"
#include "integrator.hpp"
#include <iostream>

namespace integrators {

// basic euler integrator it's the simplest but it can drift from accuracy
// easily and overshhot or undershhot aproximations it works by updating the
// position using OLD velocity (wgy it's not accurate) then updating velocity
// using acceleration
template <int DIM, typename Attributes> 
class Euler : public Integrator<DIM, Attributes> {
private:
  AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator;

public:
  Euler(AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator)
      : accelerationAccumulator(accelerationAccumulator) {}

  void step(bodies<DIM, Attributes> &bodies, double dt) override {

    accelerationAccumulator.compute(bodies);
    size_t n = bodies.localSize();
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).pos()[d] += bodies.local(i).vel()[d] * dt;
      }
    }

    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).vel()[d] += accelerationAccumulator.accel(i)[d] * dt;
      }
    }
  }
};

} // namespace integrators
