#pragma once

#include "acceleration-accumulator.hpp"
#include "forces/func.hpp"
#include <iostream>
#include <mpi.h>

template <int DIM, typename Attributes>
class MPIAccumulatorReduced : public AccelerationAccumulator<DIM, Attributes> {
private:
  const MPI_Datatype &_MPI_VEC;
  const forces::force<DIM, Attributes> &_force;
  const std::vector<double> &_massesAll;
  const std::vector<Attributes> &_attributesAll;
  int _mpiSize;
  int _mpiRank;
  std::vector<Vec<DIM>> _accelerations;
  Vec<DIM> *tmpData;
  Vec<DIM> *tmpPos;
  Vec<DIM> *tmpVel;
  Vec<DIM> *tmpAccel;
  int source;
  int dest;

  // assumes n%mpiSize==0, cyclic distribution
  static inline int global_to_local(int glb, int comm_sz) {
    return glb / comm_sz;
  }

  /* First_index(glb_part1, owner, comm_sz, n)
   * Return the first global index assigned to 'owner' that is strictly >
   * glb_part1. (Used to find which particles in the received block to pair
   * with.)
   *
   * In cyclic distribution owner has indices: owner, owner+P, owner+2P, ...
   * We find smallest glb_part2 = owner + k*P > glb_part1.
   *
   * this assumes n%mpiSize == 0
   */
  int first_index_gt(int glb_part1, int owner, int comm_sz, int n) {
    int glb = owner;
    if (glb <= glb_part1) {
      int k = (glb_part1 - owner) / comm_sz + 1;
      glb = owner + k * comm_sz;
    }
    return (glb < n) ? glb : -1;
  }

public:
  MPIAccumulatorReduced(MPI_Datatype &MPI_VEC, int localSize, int mpiSize,
                        int mpiRank, const std::vector<double> &massesAll,
                        const forces::force<DIM, Attributes> &force,
                        const std::vector<Attributes> &attributesAll)
      : _MPI_VEC(MPI_VEC), _force(force), _massesAll(massesAll),
        _attributesAll(attributesAll), _mpiSize(mpiSize), _mpiRank(mpiRank) {
    _accelerations.resize(localSize);

    tmpData = new Vec<DIM>[3 * localSize];
    tmpPos = tmpData;
    tmpVel = tmpPos + localSize;
    tmpAccel = tmpVel + localSize;

    source = (mpiRank + 1) % mpiSize;
    dest = (mpiRank - 1 + mpiSize) % mpiSize;
  }

  ~MPIAccumulatorReduced() { delete[] tmpData; }

  void compute(bodies<DIM, Attributes> &bodies) override {
    // make forces equal to zero and prepare tmpData
    for (size_t i = 0; i < bodies.localSize(); i++) {
      _accelerations[i] = 0.;
      tmpPos[i] = bodies.local(i).pos();
      tmpAccel[i] = 0.;
    }

    // first, compute forces due to interactions among local particles
    for (size_t locI = 0; locI < bodies.localSize(); locI++) {
      int glbI = _mpiRank + locI * _mpiSize;
      for (size_t locJ = locI + 1; locJ < bodies.localSize(); locJ++) {
        int glbJ = _mpiRank + locJ * _mpiSize;
        Vec<DIM> f =
            _force(bodyCopy(bodies.local(locI).pos(), bodies.local(locI).vel(),
                            _massesAll[glbI], _attributesAll[glbI]),
                   bodyCopy(bodies.local(locJ).pos(), bodies.local(locJ).vel(),
                            _massesAll[glbJ], _attributesAll[glbJ]));

        _accelerations[locI] += f;
        tmpAccel[locJ] -= f;
      }
    }

    // ring pass: mpiSize-1 phases
    for (int phase = 1; phase < _mpiSize; phase++) {
      // send tmpData (tmp Positions and tmp Forces) to dest, and receive into
      // the same buffer from the source
      MPI_Sendrecv_replace(tmpData, 3 * bodies.localSize(), _MPI_VEC, dest, 99,
                           source, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      int owner = (_mpiRank + phase) % _mpiSize;

      // compute interactions between local particles and the received ones
      for (size_t locI = 0, glbI = _mpiRank; locI < bodies.localSize();
           locI++, glbI += _mpiSize) {
        // find first global index > glbI that belongs to owner
        int glbJ = first_index_gt(glbI, owner, _mpiSize, _massesAll.size());

        while (glbJ != -1 && glbJ < glbI + _mpiSize &&
               glbJ < _massesAll.size()) {
          int locJ = global_to_local(glbJ, _mpiSize);

          // glbI < glbJ by construction
          Vec<DIM> f = _force(bodyCopy(bodies.local(locI).pos(),
                                       bodies.local(locI).vel(),
                                       _massesAll[glbI], _attributesAll[glbI]),
                              bodyCopy(tmpPos[locJ], tmpVel[locJ],
                                       _massesAll[glbJ], _attributesAll[glbJ]));

          _accelerations[locI] += f;
          tmpAccel[locJ] -= f;

          glbJ += _mpiSize;
          if (glbJ >= _massesAll.size())
            break;
        }
      }
    }

    // one last send/recv of tmp arrays
    MPI_Sendrecv_replace(tmpData, 3 * bodies.localSize(), _MPI_VEC, dest, 98,
                         source, 98, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // add the tmpForces that were just received
    for (size_t i = 0; i < bodies.localSize(); i++) {
      // FIXME: For some reason both _accelerations and tmpAccel contain opposite values.
      // Reproduce using three particles on three processes
      // Verify using:
      // std::cout << _accelerations[i][0] << " " << _accelerations[i][1] << " " << _accelerations[i][2] << " " 
      //    << tmpAccel[i][0] << " " << tmpAccel[i][1] << " " << tmpAccel[i][2] << "\n";

      _accelerations[i] += tmpAccel[i];
      for (int j = 0; j < DIM; j++)
        _accelerations[i][j] /= _massesAll[_mpiRank + i * _mpiSize];
    }
  }

  const Vec<DIM> &accel(int localIndex) const override {
    return _accelerations[localIndex];
  }
};
