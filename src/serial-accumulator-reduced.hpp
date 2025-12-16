#pragma once

#include "acceleration-accumulator.hpp"
#include "forces/func.hpp"

template <int DIM> class SerialAccumulatorReduced : public AccelerationAccumulator<DIM> {
private:
  const forces::force<DIM> &_force;
  std::vector<Vec<DIM>> _accelerations;

public:
  SerialAccumulatorReduced(int size, const forces::force<DIM> &force)
      : _force(force) {
    _accelerations.resize(size);
  }

  void compute(bodies<DIM> &bodies) override {

    size_t n = bodies.localSize();

    // azzera forze
    for (size_t i = 0; i < n; ++i)
        for (int d = 0; d < DIM; ++d)
            _accelerations[i][d] = 0.0;

    // solo coppie i < j
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            Vec<DIM> force = _force(bodies.local(i), bodies.local(j));
            _accelerations[i] += force;
            _accelerations[j] -= force;
        }
    }

    for (size_t i = 0; i < n; i++) {
      for (int j = 0; j < DIM; j++) {
        _accelerations[i][j] /= bodies.local(i).mass();
      }
    }
  }

  const Vec<DIM> &accel(int localIndex) const override {
    return _accelerations[localIndex];
  }
};
