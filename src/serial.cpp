#include "body.hpp"
#include "forces/gravity.hpp"
#include "integrators/euler.hpp"
#include "integrators/rk4.hpp"
#include "integrators/sympletic.hpp"
#include "integrators/verlet.hpp"
#ifdef REDUCED
#include "serial-accumulator-reduced.hpp"
#else
#include "serial-accumulator.hpp"
#endif
#include "utils.hpp"
#include <cmath>
#include <iostream>

constexpr int DIM = 2;

int main() {
  int n_steps;
  double dt;

  bodies<DIM, EmptyAttributes> bodies;
  utils::readFromStream(std::cin, n_steps, dt, bodies);

  forces::gravity<DIM> force;
#ifdef REDUCED
  SerialAccumulatorReduced<DIM, EmptyAttributes> accumulator(bodies.localSize(), force);
#else
  SerialAccumulator accumulator(bodies.localSize(), force);
#endif
  // integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
  integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

  // time loop
  for (int step = 0; step < n_steps; ++step) {
    integrator.step(bodies, dt);
  }

  utils::saveToStream(std::cout, n_steps, dt, bodies);

  return 0;
}
