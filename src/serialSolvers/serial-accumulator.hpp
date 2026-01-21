#pragma once

#include "../acceleration-accumulator.hpp"
#include "../forces/func.hpp"
#include "../omp_utils.hpp"

template <int DIM, typename Attributes> 
class SerialAccumulator : public AccelerationAccumulator<DIM, Attributes> {
private:
  const forces::force<DIM, Attributes> &_force;
  std::vector<Vec<DIM>> _accelerations;

public:
  SerialAccumulator(int size, const forces::force<DIM, Attributes> &force)
      : _force(force) {
    _accelerations.resize(size);
  }

  void compute(Bodies<DIM, Attributes> &bodies) override {

    size_t n = bodies.localSize();

    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i)
      _accelerations[i] = 0.0;

    OMP_STATIC_LOOP
    for (size_t q = 0; q < n; ++q) {
        for (size_t k = 0; k < n; ++k) {
            if (q == k) continue;
            _accelerations[q] += _force(bodies.local(q), bodies.local(k));
        }
    }

    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; i++) {
      _accelerations[i] /= bodies.local(i).mass();
    }
  }

  const Vec<DIM> &accel(int localIndex) const override {
    return _accelerations[localIndex];
  }
};
