#pragma once

#include "../acceleration-accumulator.hpp"
#include "../forces/func.hpp"
#include "../omp_utils.hpp"

template <int DIM, typename Attributes> 
class SerialAccumulatorReduced : public AccelerationAccumulator<DIM, Attributes> {
private:
  const forces::force<DIM, Attributes> &_force;
  std::vector<Vec<DIM>> _accelerations;

public:
  SerialAccumulatorReduced(int size, const forces::force<DIM, Attributes> &force)
      : _force(force) {
    _accelerations.resize(size);
  }

  void compute(Bodies<DIM, Attributes> &bodies) override {

    size_t n = bodies.localSize();

    OMP_STATIC_LOOP
    for (size_t i = 0; i < n; ++i)
      _accelerations[i] = 0.0;

    OMP_DYNAMIC_LOOP
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            Vec<DIM> force = _force(bodies.local(i), bodies.local(j));
            _accelerations[i] += force;
            _accelerations[j] -= force;
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
