#include "body.hpp"
#include "forces/gravity.hpp"
#ifdef REDUCED
#include "serial-accumulator-reduced.hpp"
#else
#include "serial-accumulator.hpp"
#endif
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

constexpr int DIM = 3;
constexpr double EPS = 1e-3;

double critical_timestep(const Vec<DIM> &vel, const Vec<DIM> &accel,
                         double dt_max) {
  double acc_norm_sq = 0.0;
  double vel_norm_sq = 0.0;
  for (int d = 0; d < DIM; ++d) {
    acc_norm_sq += accel[d] * accel[d];
    vel_norm_sq += vel[d] * vel[d];
  }

  double acc_norm = std::sqrt(acc_norm_sq);
  double vel_norm = std::sqrt(vel_norm_sq);

  double acc_limit = acc_norm > 0.0 ? std::sqrt(EPS / acc_norm) : dt_max;
  double vel_limit = vel_norm > 0.0 ? EPS / vel_norm : dt_max;

  double dt_crit = std::min({dt_max, acc_limit, vel_limit});
  const double min_dt = dt_max * 1e-6;
  return std::max(dt_crit, min_dt);
}

template <typename Accumulator>
double update_timesteps(const bodies<DIM, EmptyAttributes> &bodies,
                        const Accumulator &accumulator, double dt_max,
                        std::vector<double> &particle_dt) {
  const double min_dt = dt_max * 1e-6;
  double next_dt = dt_max;

  for (size_t i = 0; i < bodies.localSize(); ++i) {
    double dt_crit =
        critical_timestep(bodies.local(i).vel(), accumulator.accel(i), dt_max);
    double candidate = particle_dt[i];

    if (candidate > dt_crit) {
      while (candidate > dt_crit) {
        candidate *= 0.5; // refine until stable
      }
    } else if (candidate < 0.5 * dt_crit) {
      candidate = std::min(candidate * 2.0, dt_max); // coarsen
    }

    candidate = std::clamp(candidate, min_dt, dt_max);
    particle_dt[i] = candidate;
    next_dt = std::min(next_dt, candidate);
  }

  return next_dt;
}

int main() {
  int n;
  double total_time;
  double dt_max;

  bodies<DIM, EmptyAttributes> bodies;
  utils::readFromStream(std::cin, n, total_time, dt_max, bodies);

  forces::gravity<DIM> force;
#ifdef REDUCED
  SerialAccumulatorReduced<DIM, EmptyAttributes> accumulator(bodies.localSize(),
                                                             force);
#else
  SerialAccumulator<DIM, EmptyAttributes> accumulator(bodies.localSize(), force);
#endif

  std::vector<double> particle_dt(bodies.localSize(), dt_max);
  double time = 0.0;

  while (time < total_time) {
    accumulator.compute(bodies);

    double dt = update_timesteps(bodies, accumulator, dt_max, particle_dt);
    dt = std::min(dt, total_time - time);
    if (dt <= 0.0) {
      break;
    }

    for (size_t i = 0; i < bodies.localSize(); ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).pos()[d] += bodies.local(i).vel()[d] * dt;
        bodies.local(i).vel()[d] += accumulator.accel(i)[d] * dt;
      }
    }

    time += dt;
  }

  utils::saveToStream(std::cout, n, total_time, dt_max, bodies);

  return 0;
}
