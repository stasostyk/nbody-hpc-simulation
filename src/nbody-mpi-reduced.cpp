#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "body.hpp"
#include "forces/coulomb.hpp"
#include "forces/gravity.hpp"
#include "integrators/euler.hpp"
#include "integrators/rk4.hpp"
#include "integrators/sympletic.hpp"
#include "integrators/verlet.hpp"
#include "mpi-accumulator-reduced.hpp"
#include "utils.hpp"

#include "adaptive_dt.hpp"

constexpr int DIM = 3;

static MPI_Datatype MPI_VEC;
static MPI_Datatype MPI_CHARGE;

inline static void initMPIType() {
  MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
  MPI_Type_commit(&MPI_VEC);

  MPI_Type_contiguous(1, MPI_DOUBLE, &MPI_CHARGE);
  MPI_Type_commit(&MPI_CHARGE);
}

inline static void freeMPIType() {
  MPI_Type_free(&MPI_VEC);
  MPI_Type_free(&MPI_CHARGE);
}

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

void gatherAndSaveAllPositions(int mpiSize, int mpiRank, int n, int steps, double dt, const std::vector<Vec<DIM>> &localPositions,
                               const std::vector<Vec<DIM>> &localVelocities, const std::vector<double> &masses, const std::string &filetag) {
  bodies<DIM, EmptyAttributes> final;

  if (mpiRank == 0) {
    final.resize(n, n, 0);
  }

  for (int i = mpiRank, lidx = 0; i < n; i += mpiSize, lidx++) {
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
    for (int r = 1; r < mpiSize; r++) {
      for (int i = r; i < n; i += mpiSize) {
        MPI_Recv(&final.position[i], 1, MPI_VEC, r, 200 + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&final.velocity[i], 1, MPI_VEC, r, 400 + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }

    utils::saveToFile("test-MPI-reduced-" + filetag + ".out", steps, dt, final, false);
  }
}

void runMPIReduced(int argc, char **argv, const forces::force<DIM, EmptyAttributes> &force, int outputStride) {
  int mpiSize, mpiRank;
  allMPIInit(&argc, &argv, mpiSize, mpiRank);

  bodies<DIM, EmptyAttributes> bodies;
  std::vector<EmptyAttributes> attributes;

  int n;
  int steps;
  double dt;

  std::vector<double> masses;

  if (mpiRank == 0) {
    utils::readFromFile("test1.in.out", steps, dt, bodies);
    masses.resize(bodies.localSize());
    std::copy_n(bodies.mass.begin(), bodies.localSize(), masses.begin());
    n = bodies.globalSize();
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (n % mpiSize != 0) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int locN = n / mpiSize;

  attributes.resize(locN);

  if (mpiRank != 0) {
    masses.resize(n);
    bodies.resize(locN, locN, 0);
  }

  MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (mpiRank == 0) {
    for (int i = 0; i < n; i++) {
      int owner = i % mpiSize;
      int lidx = i / mpiSize;
      if (owner == 0) { bodies.position[lidx] = bodies.position[i]; bodies.velocity[lidx] = bodies.velocity[i];
      } else {
        MPI_Send(&bodies.position[i], 1, MPI_VEC, owner, 1, MPI_COMM_WORLD);
        MPI_Send(&bodies.velocity[i], 1, MPI_VEC, owner, 2, MPI_COMM_WORLD);
      }
    }
    bodies.resize(locN, locN, 0);
  } else {
    for (int i = mpiRank; i < n; i += mpiSize) {
      int lidx = i / mpiSize;
      MPI_Recv(&bodies.position[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&bodies.velocity[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
    }
  }

  MPIAccumulatorReduced<DIM, EmptyAttributes> accumulator(
      MPI_VEC, n / mpiSize, mpiSize, mpiRank, masses, force, attributes);
  // integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
  integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

  const double dt0 = dt;
  const double t_end = dt0 * static_cast<double>(steps);

  const double  dt_cap = dt0 * 4.0;
  double dt_max = dt_cap;

  const double dt_min = dt0 * 1e-4;
  const double eps_v = 5e-5;
  const double v_floor = 1e-6;
  const double a_floor = 1e-12;
  const double max_growth = 1.5;

  double time = 0.0;
  const int out_every = 1;
  double next_out = dt0 * out_every;
  int frame = 0;

  accumulator.compute(bodies);

  double dt_prev = dt0;

  double dt_local =
      compute_local_dt<DIM>(bodies, accumulator, dt_prev, dt_max, dt_min, eps_v,
                            v_floor, a_floor, max_growth);
  double dt_curr = dt_local;
  MPI_Allreduce(&dt_local, &dt_curr, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if (mpiRank == 0) {
    double target = std::clamp(dt_curr * 1.2, dt_min, dt_cap); // desired max
    if (target < dt_max) dt_max = target;                      // shrink fast
    else dt_max = std::min(dt_cap, dt_max * 1.05);             // grow slowly
  }
  MPI_Bcast(&dt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (mpiRank == 0) {
    std::cout << "time=" << time
              << " dt_curr=" << dt_curr
              << " dt_max=" << dt_max
              << std::endl;

  }

  while (time < t_end && !time_reached(time, t_end)) {
    double dt_step = dt_curr;
    dt_step = std::min(dt_step, t_end - time);
    // dt_step = std::min(dt_step, next_out - time);  // IMPORTANT: disable for adaptive dt


    integrator.step(bodies, dt_step);
    time += dt_step;

    accumulator.compute(bodies);

    dt_local =
        compute_local_dt<DIM>(bodies, accumulator, dt_step, dt_max, dt_min, eps_v, v_floor, a_floor, max_growth);
    MPI_Allreduce(&dt_local, &dt_curr, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    if (mpiRank == 0) {
      double target = std::clamp(dt_curr * 1.2, dt_min, dt_cap);
      if (target < dt_max) dt_max = target;
      else dt_max = std::min(dt_cap, dt_max * 1.05);
    }
    MPI_Bcast(&dt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mpiRank == 0) {
      std::cout << "time=" << time
                << " dt_step=" << dt_step
                << " dt_curr=" << dt_curr
                << " dt_max=" << dt_max
                << " dt_local(min)=" << dt_local
                << std::endl;

    }

    if (time_reached(time, next_out)) {
      if (frame % outputStride == 0) {
        gatherAndSaveAllPositions(
            mpiSize, mpiRank, n, steps, dt0, bodies.position, bodies.velocity,
            masses, std::to_string(frame / outputStride));
      }
      ++frame;
      next_out =
          std::min(t_end, dt0 * static_cast<double>((frame + 1) * out_every));
    }
  }

  gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt0, bodies.position, bodies.velocity, masses, "final");

  allMPIFinalize();
}

void run3Charges(int argc, char **argv, int outputStride) {
  int mpiSize, mpiRank;
  allMPIInit(&argc, &argv, mpiSize, mpiRank);

  bodies<DIM, forces::charge> bodies;

  int n;
  int steps;
  double dt;

  std::vector<double> masses;
  std::vector<forces::charge> attributes;

  if (mpiRank == 0) {
    steps = 1000;
    dt = 0.01;

    bodies.resize(3, 3, 0);
    bodies.local(0).pos()[0] = 1.0;
    bodies.local(0).pos()[1] = 1.0;
    bodies.local(0).pos()[2] = 1.0;
    bodies.local(0).vel()[0] = 0.0;
    bodies.local(0).vel()[1] = 0.0;
    bodies.local(0).vel()[2] = 0.0;
    bodies.local(0).mass() = 1.0;
    bodies.local(0).attributes().charge = 1e-4;

    bodies.local(1).pos()[0] = -1.0;
    bodies.local(1).pos()[1] = -1.0;
    bodies.local(1).pos()[2] = 1.0;
    bodies.local(1).vel()[0] = 0.0;
    bodies.local(1).vel()[1] = 0.0;
    bodies.local(1).vel()[2] = 0.0;
    bodies.local(1).mass() = 1.0;
    bodies.local(1).attributes().charge = 1e-4;

    bodies.local(2).pos()[0] = -1.0;
    bodies.local(2).pos()[1] = 1.0;
    bodies.local(2).pos()[2] = 1.0;
    bodies.local(2).vel()[0] = 0.0;
    bodies.local(2).vel()[1] = 0.0;
    bodies.local(2).vel()[2] = 0.0;
    bodies.local(2).mass() = 1.0;
    bodies.local(2).attributes().charge = 1e-4;

    masses.resize(bodies.localSize());
    std::copy_n(bodies.mass.begin(), bodies.localSize(), masses.begin());
    attributes.resize(bodies.localSize());
    std::copy_n(bodies.attributes.begin(), bodies.localSize(), attributes.begin());

    n = bodies.globalSize();
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (n % mpiSize != 0) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int locN = n / mpiSize;

  if (mpiRank != 0) {
    masses.resize(n);
    attributes.resize(n);
    bodies.resize(locN, locN, 0);
  }

  MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(attributes.data(), n, MPI_CHARGE, 0, MPI_COMM_WORLD);

  if (mpiRank == 0) {
    for (int i = 0; i < n; i++) {
      int owner = i % mpiSize;
      int lidx = i / mpiSize;
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
      int lidx = i / mpiSize;
      MPI_Recv(&bodies.position[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&bodies.velocity[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  forces::coulomb<DIM> coulomb{};
  const forces::force<DIM, forces::charge> &force_c = coulomb;

  MPIAccumulatorReduced<DIM, forces::charge> accumulator(
      MPI_VEC, n / mpiSize, mpiSize, mpiRank, masses, force_c, attributes);
  // integrators::Euler<DIM, forces::charge> integrator(accumulator);
  // integrators::Sympletic<DIM, forces::charge> integrator(accumulator);
  // integrators::Verlet<DIM, forces::charge> integrator(accumulator);
  integrators::RK4<DIM, forces::charge> integrator(accumulator);

  for (int step = 0; step < steps; step++) {
    integrator.step(bodies, dt);

    if (step % outputStride == 0) {
      gatherAndSaveAllPositions(
          mpiSize, mpiRank, n, steps, dt, bodies.position, bodies.velocity,
          masses, std::to_string(step / outputStride));
    }
  }

  gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt, bodies.position, bodies.velocity, masses, "final");

  allMPIFinalize();
}

int main(int argc, char **argv) {
  forces::gravity<DIM> gravity{};
  const forces::force<DIM, EmptyAttributes> &force = gravity;

  runMPIReduced(argc, argv, force, 20);
  return 0;
}
