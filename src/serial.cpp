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
#include "timer.hpp"
#include <cmath>
#include <iostream>

constexpr int DIM = 2;

int main() {
  int n_steps;
  double dt;

  bodies<DIM, EmptyAttributes> bodies;
  // utils::readFromStream(std::cin, n_steps, dt, bodies);
  utils::readFromFile<DIM>("test1.in.out", n_steps, dt, bodies);

  forces::gravity<DIM> force;
#ifdef REDUCED
  SerialAccumulatorReduced<DIM, EmptyAttributes> accumulator(bodies.localSize(), force);
#else
  SerialAccumulator accumulator(bodies.localSize(), force);
#endif
  integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

  Timer t;
  t.start();

  // time loop
  for (int step = 0; step < n_steps; ++step) {
    integrator.step(bodies, dt);
  }

  t.end();
  t.print();

  // utils::saveToStream(std::cout, n_steps, dt, bodies);
  utils::saveToFile<DIM>("serial.test1.lastStep.out", n_steps, dt, bodies, false);

  utils::compareOutputsSingleFile<DIM>("bhSerial.test1.lastStep.out", "serial.test1.lastStep.out");

  return 0;
}
