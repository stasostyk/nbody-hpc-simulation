#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "bodies.hpp"
#include "vec.hpp"

template <int DIM>
inline double vec_norm(const Vec<DIM> &v) {
  double s = 0.0;
  for (int d = 0; d < DIM; ++d) s += v[d] * v[d];
  return std::sqrt(s);
}

template <int DIM, class Bodies, class Acc>
double compute_local_dt_raw(const Bodies &b, const Acc &acc,
                            double eps_v, double v_floor, double a_floor)
{
  double local_min = 1e100;
  for (int i = 0; i < b.localSize(); ++i) {
    const auto v = b.velocity[i + b.localOffset()];
    const auto a = acc.accel(i);
    const double vnorm = std::sqrt(v.dot(v));
    const double anorm = std::sqrt(a.dot(a));
    const double dt_i = eps_v * (std::max(vnorm, v_floor) / std::max(anorm, a_floor));
    local_min = std::min(local_min, dt_i);
  }
  return local_min;
}


template <int DIM, typename Attributes, typename Accumulator>
inline double compute_local_dt(const bodies<DIM, Attributes> &b, const Accumulator &acc,
                               double dt_prev, double dt_max, double dt_min,
                               double eps_v, double v_floor, double a_floor,
                               double max_growth) {
  double local_min = std::numeric_limits<double>::infinity();
  const size_t n = b.localSize();
  for (size_t i = 0; i < n; ++i) {
    const double v = vec_norm<DIM>(b.local(i).vel());
    const double a = vec_norm<DIM>(acc.accel((int)i));
    const double denom = std::max(a, a_floor);
    const double numer = std::max(v, v_floor);
    const double dt_i = eps_v * (numer / denom);
    local_min = std::min(local_min, dt_i);
  }

  double dt = local_min;
  if (!std::isfinite(dt)) dt = dt_max;
  dt = std::min(dt, dt_max);
  if (dt_prev > 0.0) dt = std::min(dt, dt_prev * max_growth);
  dt = std::max(dt, dt_min);
  return dt;
}

inline bool time_reached(double t, double target) {
  const double tol = 1e-12 * std::max(1.0, std::abs(target));
  return std::abs(t - target) <= tol;
}
