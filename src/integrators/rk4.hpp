#pragma once

#include "../acceleration-accumulator.hpp"
#include "../body.hpp"
#include "integrator.hpp"
#include <algorithm>

namespace integrators {

// Runge-Kutta 4 integrator (RK4) is much more accurate than Euler and Verlet
// it estimates the slope (derivative) four times per timestep
// and then takes a weighted average of those slopes
// this gives very good accuracy for short and medium time simulations
// but it's not symplectic, so energy can still slowly drift in long runs
template <int DIM, typename Attributes> 
class RK4 : public Integrator<DIM, Attributes> {
private:
  AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator;

public:
  RK4(AccelerationAccumulator<DIM, Attributes>& accelerationAccumulator)
      : accelerationAccumulator(accelerationAccumulator) {}

public:
  void step(Bodies<DIM, Attributes> &bodies, double dt) override {

    size_t n = bodies.localSize();
    // k?_s are the slopes for position (ds/dt = v)
    // k?_v are the slopes for velocity (dv/dt = a)
    std::vector<Vec<DIM>> k1_s(n), k2_s(n), k3_s(n), k4_s(n);
    std::vector<Vec<DIM>> k1_v(n), k2_v(n), k3_v(n), k4_v(n);
    // Temporary copy for intermediate steps
    std::vector<Vec<DIM>> posOriginal(
        bodies.position.begin() + bodies.localOffset(),
        bodies.position.begin() + bodies.localOffset() + n);
    std::vector<Vec<DIM>> velOriginal(
        bodies.velocity.begin() + bodies.localOffset(),
        bodies.velocity.begin() + bodies.localOffset() + n);

    // Step 1: Compute k1
    accelerationAccumulator.compute(bodies);
    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      k1_s[i] = bodies.local(i).vel();
      k1_v[i] = accelerationAccumulator.accel(i);
    }

    // Step 2: Compute k2
    // we move the system forward by dt/2 using k1 slopes
    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      bodies.local(i).pos() = posOriginal[i] + (0.5 * dt) * k1_s[i];
      bodies.local(i).vel() = velOriginal[i] + (0.5 * dt) * k1_v[i];
    }

    accelerationAccumulator.compute(bodies); // a at s^(2)

    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      k2_s[i] = bodies.local(i).vel();
      k2_v[i] = accelerationAccumulator.accel(i);
    }

    // Step 3: compute k3, again at dt/2 but now using k2 slopes
    // this improves the midpoint estimate
    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      bodies.local(i).pos() = posOriginal[i] + (0.5 * dt) * k2_s[i];
      bodies.local(i).vel() = velOriginal[i] + (0.5 * dt) * k2_v[i];
    }

    accelerationAccumulator.compute(bodies);

    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      k3_s[i] = bodies.local(i).vel();
      k3_v[i] = accelerationAccumulator.accel(i);
    }

    // Step 4: compute k4 using a full timestep with k3
    // this gives the slope at the end of the timestep
    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      bodies.local(i).pos() = posOriginal[i] + dt * k3_s[i];
      bodies.local(i).vel() = velOriginal[i] + dt * k3_v[i];
    }

    accelerationAccumulator.compute(bodies);

    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      k4_s[i] = bodies.local(i).vel();
      k4_v[i] = accelerationAccumulator.accel(i);
    }

    std::copy_n(posOriginal.begin(), n, bodies.position.begin() + bodies.localOffset());
    std::copy_n(velOriginal.begin(), n, bodies.velocity.begin() + bodies.localOffset());

    // Final and actual update
    // we combine all four slopes using RK4 weights:
    // k1 + 2*k2 + 2*k3 + k4, divided by 6
    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i) {
      bodies.local(i).vel() += (dt / 6.0) * (k1_v[i] + 2.0 * k2_v[i] +
                                          2.0 * k3_v[i] + k4_v[i]);
      bodies.local(i).pos() += (dt / 6.0) * (k1_s[i] + 2.0 * k2_s[i] +
                                          2.0 * k3_s[i] + k4_s[i]);
    }
  }
};

} // namespace integrators
