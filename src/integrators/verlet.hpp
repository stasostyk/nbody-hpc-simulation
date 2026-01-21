#pragma once

#include "../acceleration-accumulator.hpp"
#include "integrator.hpp"

namespace integrators {

// velocity verlet integrator is even more accurate
// it updates velocity in half steps and then use that half step v to compute
// position update recomputes acceleration using updated positions this way it
// takes into account the change in acceleration due to position change
template <int DIM, typename Attributes> 
class Verlet : public Integrator<DIM, Attributes> {
private:
  AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator;

public:
  Verlet(AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator)
      : accelerationAccumulator(accelerationAccumulator) {}

  void step(bodies<DIM, Attributes> &bodies, double dt) override {

    accelerationAccumulator.compute(bodies);
    size_t n = bodies.localSize();
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        // we update velocity only by half step
        bodies.local(i).vel()[d] += accelerationAccumulator.accel(i)[d] * 0.5 * dt;
        // then we update position using this half step velocity
        bodies.local(i).pos()[d] += bodies.local(i).vel()[d] * dt;
      }
    }
    // recompute acceleration with new positions
    accelerationAccumulator.compute(bodies);
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        // complete velocity update with the new acceleration
        bodies.local(i).vel()[d] += accelerationAccumulator.accel(i)[d] * 0.5 * dt;
      }
    }
  }
};

} // namespace integrators
