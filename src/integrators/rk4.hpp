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
template <int DIM> class RK4 : public Integrator<DIM> {
private:
  AccelerationAccumulator<DIM>& accelerationAccumulator;

public:
  RK4(AccelerationAccumulator<DIM>& accelerationAccumulator)
      : accelerationAccumulator(accelerationAccumulator) {}

public:
  void step(bodies<DIM> &bodies, double dt) override {

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
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        k1_s[i][d] = bodies.local(i).vel()[d];
        k1_v[i][d] = accelerationAccumulator.accel(i)[d];
      }
    }

    // Step 2: Compute k2
    // we move the system forward by dt/2 using k1 slopes
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).pos()[d] = posOriginal[i][d] + (0.5 * dt) * k1_s[i][d];
        bodies.local(i).vel()[d] = velOriginal[i][d] + (0.5 * dt) * k1_v[i][d];
      }
    }

    accelerationAccumulator.compute(bodies); // a at s^(2)

    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        k2_s[i][d] = bodies.local(i).vel()[d];
        k2_v[i][d] = accelerationAccumulator.accel(i)[d];
      }
    }

    // Step 3: compute k3, again at dt/2 but now using k2 slopes
    // this improves the midpoint estimate
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).pos()[d] = posOriginal[i][d] + (0.5 * dt) * k2_s[i][d];
        bodies.local(i).vel()[d] = velOriginal[i][d] + (0.5 * dt) * k2_v[i][d];
      }
    }

    accelerationAccumulator.compute(bodies);

    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        k3_s[i][d] = bodies.local(i).vel()[d];
        k3_v[i][d] = accelerationAccumulator.accel(i)[d];
      }
    }

    // Step 4: compute k4 using a full timestep with k3
    // this gives the slope at the end of the timestep
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).pos()[d] = posOriginal[i][d] + dt * k3_s[i][d];
        bodies.local(i).vel()[d] = velOriginal[i][d] + dt * k3_v[i][d];
      }
    }

    accelerationAccumulator.compute(bodies);

    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        k4_s[i][d] = bodies.local(i).vel()[d];
        k4_v[i][d] = accelerationAccumulator.accel(i)[d];
      }
    }

    std::copy_n(posOriginal.begin(), n, bodies.position.begin() + bodies.localOffset());
    std::copy_n(velOriginal.begin(), n, bodies.velocity.begin() + bodies.localOffset());

    // Final and actual update
    // we combine all four slopes using RK4 weights:
    // k1 + 2*k2 + 2*k3 + k4, divided by 6
    for (size_t i = 0; i < n; ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).vel()[d] += (dt / 6.0) * (k1_v[i][d] + 2.0 * k2_v[i][d] +
                                            2.0 * k3_v[i][d] + k4_v[i][d]);
        bodies.local(i).pos()[d] += (dt / 6.0) * (k1_s[i][d] + 2.0 * k2_s[i][d] +
                                            2.0 * k3_s[i][d] + k4_s[i][d]);
      }
    }
  }
};

} // namespace integrators
