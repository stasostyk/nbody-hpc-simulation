#pragma once

#include "acceleration-accumulator.hpp"
#include "forces/func.hpp"
#include <mpi.h>

template <int DIM, typename Attributes> 
class MPIAccumulator : public AccelerationAccumulator<DIM, Attributes> {
private:
  const MPI_Datatype &_MPI_VEC;
  const forces::force<DIM, Attributes> &_force;
  const std::vector<int> &_counts;
  const std::vector<int> &_offsets;
  std::vector<Vec<DIM>> _accelerations;

public:
  MPIAccumulator(MPI_Datatype &MPI_VEC, int localSize,
                 const std::vector<int> &counts,
                 const std::vector<int> &offsets,
                 const forces::force<DIM, Attributes> &force)
      : _MPI_VEC(MPI_VEC), _force(force), _counts(counts), _offsets(offsets) {
    _accelerations.resize(localSize);
  }

  void compute(Bodies<DIM, Attributes> &bodies) override {
    // All tasks gather all positions (will be required for the next step)
    MPI_Allgatherv(bodies.position.data() + bodies.localOffset(),
                   bodies.localSize(), _MPI_VEC, bodies.position.data(),
                   _counts.data(), _offsets.data(), _MPI_VEC, MPI_COMM_WORLD);

    // Considering a hypothetical force that depends on velocity, it also must
    // be gathered
    MPI_Allgatherv(bodies.velocity.data() + bodies.localOffset(),
                   bodies.localSize(), _MPI_VEC, bodies.velocity.data(),
                   _counts.data(), _offsets.data(), _MPI_VEC, MPI_COMM_WORLD);

    // make forces equal to zero
    for (size_t i = 0; i < bodies.localSize(); i++) {
      for (int j = 0; j < DIM; j++) {
        _accelerations[i][j] = 0.;
      }
    }

    // compute all forces
    for (size_t i = 0; i < bodies.localSize(); i++) {
      size_t q = bodies.localOffset() + i;
      for (size_t k = 0; k < bodies.globalSize(); k++) {
        if (q == k)
          continue;
        _accelerations[i] += _force(bodies.local(i), bodies.global(k));
      }
    }

    for (size_t i = 0; i < bodies.localSize(); i++) {
      for (int j = 0; j < DIM; j++) {
        _accelerations[i][j] /= bodies.local(i).mass();
      }
    }
  }

  const Vec<DIM> &accel(int localIndex) const override {
    return _accelerations[localIndex];
  }
};
