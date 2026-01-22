#pragma once

#include "BHTree.hpp"
#include "../acceleration-accumulator.hpp"
#include "../forces/func.hpp"

#if USE_OPENMP
    #include <omp.h>
#endif

template <int DIM, typename Attributes> 
class BarnesHutAccumulator : public AccelerationAccumulator<DIM, Attributes> {
private:
  const int n;
  const double theta;
  const forces::force<DIM, Attributes> &_force;
  std::vector<Vec<DIM>> _accelerations;
  const Box<DIM> universeBounds;

public:
  BarnesHutAccumulator(int _n, double _theta, const Box<DIM> &_universeBounds,
                       const forces::force<DIM, Attributes> &force)
      : n(_n), theta(_theta), _force(force), universeBounds(_universeBounds) {
    _accelerations.resize(n);
  }

  void compute(bodies<DIM, Attributes> &bodiesPassed) override {
    std::vector<Body<DIM>> bodies(n);
    for (int i = 0; i < n; i++) {
        bodies[i].bodyId = i;
        bodies[i].acceleration = 0.;
        bodies[i].force = 0.;
        bodies[i].mass = bodiesPassed.global(i).mass();
        bodies[i].position = bodiesPassed.global(i).pos();
        bodies[i].velocity = bodiesPassed.global(i).vel();
    }

    // Construct tree
    BHTree<DIM, Attributes> bhTree(universeBounds, bodies);

    // Compute all forces using BH Tree.
    #if USE_OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
    #endif
    for (int i = 0; i < n; i++) {
        _accelerations[i] = bhTree.calculateForce(bodies[i], theta, _force);
    }

    for (size_t i = 0; i < (size_t)n; i++) {
      for (int j = 0; j < DIM; j++) {
        _accelerations[i][j] /= bodiesPassed.global(i).mass();
      }
    }
  }

  const Vec<DIM> &accel(int localIndex) const override {
    return _accelerations[localIndex];
  }
};
