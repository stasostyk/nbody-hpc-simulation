#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "body.hpp"

namespace timestep {

constexpr double EPS = 1e-3;

template <int DIM, typename Attributes>
inline double critical_timestep(const body<DIM, Attributes> &b,
                                const Vec<DIM> &accel,
                                double dt_max)
{
  double acc_norm_sq = 0.0;
  double vel_norm_sq = 0.0;

  for (int d = 0; d < DIM; ++d) {
    acc_norm_sq += accel[d] * accel[d];
    vel_norm_sq += b.vel()[d] * b.vel()[d];
  }

  const double acc_norm = std::sqrt(acc_norm_sq);
  const double vel_norm = std::sqrt(vel_norm_sq);

  const double acc_limit = (acc_norm > 0.0) ? std::sqrt(EPS / acc_norm) : dt_max;
  const double vel_limit = (vel_norm > 0.0) ? (EPS / vel_norm) : dt_max;

  double dt_crit = std::min({dt_max, acc_limit, vel_limit});

  const double dt_min = dt_max * 1e-6;
  return std::max(dt_crit, dt_min);
}

template <int DIM, typename Attributes, typename AccT>
inline double update_timesteps(const bodies<DIM, Attributes> &bodies,
                               const AccT &acc,
                               double dt_max,
                               std::vector<double> &particle_dt)
{
  const double dt_min = dt_max * 1e-6;

  if (particle_dt.size() != bodies.localSize())
    particle_dt.assign(bodies.localSize(), dt_max);

  double next_dt = dt_max;

  for (size_t i = 0; i < bodies.localSize(); ++i) {
    const double dt_crit =
        critical_timestep<DIM, Attributes>(bodies.local((int)i),
                                           acc.accel((int)i),
                                           dt_max);

    double candidate = particle_dt[i];

    if (candidate > dt_crit) {
      while (candidate > dt_crit)
        candidate *= 0.5;
    } else if (candidate < 0.5 * dt_crit) {
      candidate = std::min(candidate * 2.0, dt_max);
    }

    candidate = std::clamp(candidate, dt_min, dt_max);
    particle_dt[i] = candidate;

    next_dt = std::min(next_dt, candidate);
  }

  return next_dt;
}

}
