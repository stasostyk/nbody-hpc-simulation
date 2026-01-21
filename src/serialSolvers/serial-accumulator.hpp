#pragma once

#include "../acceleration-accumulator.hpp"
#include "../forces/func.hpp"

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

    // azzera forze
    for (size_t i = 0; i < n; ++i)
        for (int d = 0; d < DIM; ++d)
            _accelerations[i][d] = 0.0;

    // tutte le coppie q != k
    for (size_t q = 0; q < n; ++q) {
        for (size_t k = 0; k < n; ++k) {
            if (q == k) continue;
            _accelerations[q] += _force(bodies.local(q), bodies.local(k));
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
