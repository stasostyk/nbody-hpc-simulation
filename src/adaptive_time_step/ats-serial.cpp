#include "../body.hpp"
#include "../forces/gravity.hpp"
#include "../integrators/euler.hpp"
#include "../integrators/rk4.hpp"
#include "../integrators/sympletic.hpp"
#include "../integrators/verlet.hpp"
#ifdef REDUCED
#include "../serialSolvers/serial-accumulator-reduced.hpp"
#else
#include "../serialSolvers/serial-accumulator.hpp"
#endif
#include "../utils.hpp"
#include "../timer.hpp"
#include "adaptive_dt.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

constexpr int DIM = 3;

int main() {
  int n_steps;
  double dt;

  Bodies<DIM, EmptyAttributes> bodies;
  utils::readFromFile<DIM>("test1.in.out", n_steps, dt, bodies);

  forces::gravity<DIM> force;
#ifdef REDUCED
  SerialAccumulatorReduced<DIM, EmptyAttributes> accumulator(bodies.localSize(), force);
#else
  SerialAccumulator accumulator(bodies.localSize(), force);
#endif

  integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);

  Timer t;
  t.start();

  const double dt0 = dt;
  const double t_end = dt0 * static_cast<double>(n_steps);

  const double dt_max = dt0 * 4;
  const double dt_min = dt0 * 1e-4;
  const double eps_v = 5e-5;
  const double v_floor = 1e-6;
  const double a_floor = 1e-12;
  const double max_growth = 1.5;

  double time = 0.0;
  const int out_every = 1;
  double next_out = dt0 * out_every;
  int frame = 0;

  accumulator.compute(bodies);

  double dt_step = dt0;

  const double dt_raw0 =
    compute_local_dt<DIM>(bodies, accumulator, dt_step, 1e100, dt_min, eps_v, v_floor, a_floor, 1e100);

  double dt_curr =
    compute_local_dt<DIM>(bodies, accumulator, dt_step, dt_max, dt_min, eps_v, v_floor, a_floor, max_growth);

  std::cout << "time=" << time << " dt_raw=" << dt_raw0 << " dt_curr=" << dt_curr << " dt_max=" << dt_max << std::endl;

  while (time < t_end && !time_reached(time, t_end)) {
  dt_step = dt_curr;

  const double to_end = std::max(0.0, t_end - time);
  //const double to_out = std::max(0.0, next_out - time);

  dt_step = std::min(dt_step, to_end);
  //dt_step = std::min(dt_step, to_out);

  // Optional: guard
  if (dt_step <= 0.0) {
    if (time_reached(time, next_out)) {
      ++frame;
      next_out = std::min(t_end, dt0 * static_cast<double>((frame + 1) * out_every));
      continue;
    }
    break;
  }

  integrator.step(bodies, dt_step);
  time += dt_step;

  accumulator.compute(bodies);

  const double dt_raw =
    compute_local_dt<DIM>(bodies, accumulator, dt_step, 1e100, dt_min, eps_v, v_floor, a_floor, 1e100);

  dt_curr =
    compute_local_dt<DIM>(bodies, accumulator, dt_step, dt_max, dt_min, eps_v, v_floor, a_floor, max_growth);

  std::cout << "time=" << time << " dt_step=" << dt_step
            << " dt_curr=" << dt_curr << " dt_raw=" << dt_raw << std::endl;

  if (time_reached(time, next_out)) {
    ++frame;
    next_out = std::min(t_end, dt0 * static_cast<double>((frame + 1) * out_every));
  }
}

  t.end();
  t.print();

  utils::saveToFile<DIM>("serial.test1.lastStep.out", n_steps, dt0, bodies, false);
  //utils::compareOutputsSingleFile<DIM>("bhSerial.test1.lastStep.out", "serial.test1.lastStep.out");
  return 0;
}
