#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "acceleration-accumulator.hpp"
#include "body.hpp"
#include "forces/gravity.hpp"
#include "mpi-accumulator.hpp"
#include "utils.hpp"

constexpr int DIM = 3;
constexpr double EPS = 1e-3;

static MPI_Datatype MPI_VEC;

inline static void initMPIType() {
  MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
  MPI_Type_commit(&MPI_VEC);
}

inline static void freeMPIType() { MPI_Type_free(&MPI_VEC); }

// utility: block counts/displacements for an array of N elements across comm_sz
void make_counts_displs(int N, int comm_sz, std::vector<int> &counts,
                        std::vector<int> &displs) {
  counts.resize(comm_sz);
  displs.resize(comm_sz);
  int base = N / comm_sz;
  int rem = N % comm_sz;
  int off = 0;
  for (int i = 0; i < comm_sz; ++i) {
    counts[i] = base + (i < rem ? 1 : 0);
    displs[i] = off;
    off += counts[i];
  }
}

inline void allMPIInit(int *argc, char ***argv, int &mpiSize, int &mpiRank) {
  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  initMPIType(); // to be able to use MPI_VEC for Vecs
}

inline void allMPIFinalize() {
  freeMPIType();
  MPI_Finalize();
}

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
        candidate *= 0.5;
      }
    } else if (candidate < 0.5 * dt_crit) {
      candidate = std::min(candidate * 2.0, dt_max);
    }

    candidate = std::clamp(candidate, min_dt, dt_max);
    particle_dt[i] = candidate;
    next_dt = std::min(next_dt, candidate);
  }

  return next_dt;
}

int main(int argc, char **argv) {
  constexpr int outputStride = 1;

  const forces::force<DIM, EmptyAttributes> &force = forces::gravity<DIM>();
  bodies<DIM, EmptyAttributes> bodies;

  int mpiSize, mpiRank;
  allMPIInit(&argc, &argv, mpiSize, mpiRank);

  int n;
  double total_time;
  double dt_max;

  if (mpiRank == 0) {
    // Root process reads input, and will broadcast the data.
    utils::readFromFile("test1.in.out", n, total_time, dt_max, bodies);
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&total_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // TODO make this computation only in root process (mpiRank 0)
  std::vector<int> counts, displs;
  make_counts_displs(n, mpiSize, counts, displs);
  int locN = counts[mpiRank];
  int locOffset = displs[mpiRank];

  bodies.resize(n, locN, locOffset);

  MPI_Bcast(bodies.position.data(), n, MPI_VEC, 0, MPI_COMM_WORLD);
  MPI_Bcast(bodies.velocity.data(), n, MPI_VEC, 0, MPI_COMM_WORLD);
  MPI_Bcast(bodies.mass.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  MPIAccumulator<DIM, EmptyAttributes> accumulator(MPI_VEC, locN, counts,
                                                    displs, force);
  std::vector<double> particle_dt(locN, dt_max);

  double time = 0.0;
  int step = 0;

  while (time < total_time) {
    accumulator.compute(bodies);

    double local_next_dt =
        update_timesteps(bodies, accumulator, dt_max, particle_dt);
    double next_dt = dt_max;
    MPI_Allreduce(&local_next_dt, &next_dt, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);

    double dt = std::min(next_dt, total_time - time);
    if (dt <= 0.0) {
      break;
    }

    for (size_t i = 0; i < bodies.localSize(); ++i) {
      for (int d = 0; d < DIM; ++d) {
        bodies.local(i).pos()[d] += bodies.local(i).vel()[d] * dt;
        bodies.local(i).vel()[d] += accumulator.accel(i)[d] * dt;
      }
    }

    if (step % outputStride == 0) {
      MPI_Allgatherv(bodies.position.data() + bodies.localOffset(),
                     bodies.localSize(), MPI_VEC, bodies.position.data(),
                     counts.data(), displs.data(), MPI_VEC, MPI_COMM_WORLD);
      MPI_Allgatherv(bodies.velocity.data() + bodies.localOffset(),
                     bodies.localSize(), MPI_VEC, bodies.velocity.data(),
                     counts.data(), displs.data(), MPI_VEC, MPI_COMM_WORLD);

      if (mpiRank == 0) {
        // save for testing
        utils::saveToFile("test1-MPI." + std::to_string(step) + ".out", n,
                          total_time, dt_max, bodies, false);
      }
    }

    time += dt;
    ++step;
  }

  allMPIFinalize();
}
