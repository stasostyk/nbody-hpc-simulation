#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <string>

#include "body.hpp"
#include "forces/gravity.hpp"
#include "integrators/rk4.hpp"
#include "mpi-accumulator-reduced.hpp"
#include "utils.hpp"
#include "adaptive_dt.hpp"

constexpr int DIM = 3;

static MPI_Datatype MPI_VEC;

inline static void initMPIType() {
  MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
  MPI_Type_commit(&MPI_VEC);
}

inline static void freeMPIType() { MPI_Type_free(&MPI_VEC); }

inline void allMPIInit(int *argc, char ***argv, int &mpiSize, int &mpiRank) {
  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  initMPIType();
}

inline void allMPIFinalize() {
  freeMPIType();
  MPI_Finalize();
}

static void gatherAndSaveAllPositions(int mpiSize, int mpiRank, int n, int steps,
                                      double dt,
                                      const std::vector<Vec<DIM>> &localPositions,
                                      const std::vector<Vec<DIM>> &localVelocities,
                                      const std::vector<double> &masses,
                                      const std::string &filetag) {
  bodies<DIM, EmptyAttributes> final;
  if (mpiRank == 0) {
    final.resize(n, n, 0);
  }

  for (int i = mpiRank, lidx = 0; i < n; i += mpiSize, ++lidx) {
    if (mpiRank == 0) {
      final.position[i] = localPositions[lidx];
      final.velocity[i] = localVelocities[lidx];
    } else {
      MPI_Send(&localPositions[lidx], 1, MPI_VEC, 0, 200 + i, MPI_COMM_WORLD);
      MPI_Send(&localVelocities[lidx], 1, MPI_VEC, 0, 400 + i, MPI_COMM_WORLD);
    }
  }

  if (mpiRank == 0) {
    std::copy_n(masses.begin(), n, final.mass.begin());
    for (int r = 1; r < mpiSize; ++r) {
      for (int i = r; i < n; i += mpiSize) {
        MPI_Recv(&final.position[i], 1, MPI_VEC, r, 200 + i, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(&final.velocity[i], 1, MPI_VEC, r, 400 + i, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }
    }

    utils::saveToFile("test-MPI-reduced-" + filetag + ".out", steps, dt, final, false);
  }
}

static void runMPIReduced(int argc, char **argv,
                          const forces::force<DIM, EmptyAttributes> &force,
                          int outputStride) {
  int mpiSize, mpiRank;
  allMPIInit(&argc, &argv, mpiSize, mpiRank);

  bodies<DIM, EmptyAttributes> bodies;

  int n = 0;
  int steps = 0;
  double dt = 0.0;

  std::vector<double> masses;
  std::vector<EmptyAttributes> attributes;

  if (mpiRank == 0) {
    utils::readFromFile("test1.in.out", steps, dt, bodies);
    n = bodies.globalSize();
    masses.resize(bodies.localSize());
    std::copy_n(bodies.mass.begin(), bodies.localSize(), masses.begin());
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  attributes.resize(n);

  if (n % mpiSize != 0) {
    if (mpiRank == 0) {
      std::cout << "n should be divisible by mpiSize\n";
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  const int locN = n / mpiSize;

  if (mpiRank != 0) {
    masses.resize(n);
    bodies.resize(locN, locN, 0);
  }

  MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // cyclic scatter: i -> owner (i % mpiSize), local index = i / mpiSize
  if (mpiRank == 0) {
    for (int i = 0; i < n; ++i) {
      const int owner = i % mpiSize;
      const int lidx  = i / mpiSize;
      if (owner == 0) {
        bodies.position[lidx] = bodies.position[i];
        bodies.velocity[lidx] = bodies.velocity[i];
      } else {
        MPI_Send(&bodies.position[i], 1, MPI_VEC, owner, 1, MPI_COMM_WORLD);
        MPI_Send(&bodies.velocity[i], 1, MPI_VEC, owner, 2, MPI_COMM_WORLD);
      }
    }
    bodies.resize(locN, locN, 0);
  } else {
    for (int i = mpiRank; i < n; i += mpiSize) {
      const int lidx = i / mpiSize;
      MPI_Recv(&bodies.position[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&bodies.velocity[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }

  MPIAccumulatorReduced<DIM, EmptyAttributes> accumulator(
      MPI_VEC, locN, mpiSize, mpiRank, masses, force, attributes);

  integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

  const double dt0   = dt;
  const double t_end = dt0 * static_cast<double>(steps);

  const double dt_max     = dt0 * 4.0;
  const double dt_min     = dt0 * 1e-6;
  const double eps_v      = 5e-7;
  const double v_floor    = 1e-6;
  const double a_floor    = 1e-12;
  const double max_growth = 1.5;

  double time = 0.0;
  const int out_every = 10;
  double next_out = dt0 * out_every;
  int frame = 0;

  accumulator.compute(bodies);

  double dt_local = compute_local_dt<DIM>(bodies, accumulator, dt0,
                                         dt_max, dt_min, eps_v, v_floor, a_floor, max_growth);
  double dt_curr = dt_local;
  MPI_Allreduce(&dt_local, &dt_curr, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  while (time < t_end && !time_reached(time, t_end)) {
    double dt_step = dt_curr;
    dt_step = std::min(dt_step, t_end - time);
    dt_step = std::min(dt_step, next_out - time);

    integrator.step(bodies, dt_step);
    time += dt_step;

    accumulator.compute(bodies);

    dt_local = compute_local_dt<DIM>(bodies, accumulator, dt_step,
                                     dt_max, dt_min, eps_v, v_floor, a_floor, max_growth);
    MPI_Allreduce(&dt_local, &dt_curr, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    if (time_reached(time, next_out)) {
      if (frame % outputStride == 0) {
        gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt0,
                                  bodies.position, bodies.velocity, masses,
                                  std::to_string(frame / outputStride));
      }
      ++frame;
      next_out = std::min(t_end, dt0 * static_cast<double>((frame + 1) * out_every));
    }
  }

  gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt0,
                            bodies.position, bodies.velocity, masses, "final");

  allMPIFinalize();
}

int main(int argc, char **argv) {
  forces::gravity<DIM> gravity{};
  const forces::force<DIM, EmptyAttributes> &force = gravity;
  runMPIReduced(argc, argv, force, 20);
  return 0;
}
