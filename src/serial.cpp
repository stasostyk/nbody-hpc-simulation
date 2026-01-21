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
#include "adaptive_timestep.hpp"
#include <vector>
#include <algorithm>


constexpr int DIM = 3;

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

  const double dt_max = dt;
  const double T_end  = n_steps * dt_max;

  double time = 0.0;
  std::vector<double> particle_dt(bodies.localSize(), dt_max);

  while (time + 1e-15 < T_end) {
    accumulator.compute(bodies);

    const double dt_local =
        timestep::update_timesteps<DIM, EmptyAttributes>(
            bodies, accumulator, dt_max, particle_dt);

    const double dt_step = std::min(dt_local, T_end - time);

    integrator.step(bodies, dt_step);
    time += dt_step;
  }


  utils::saveToStream(std::cout, n_steps, dt, bodies);
  return 0;
}
