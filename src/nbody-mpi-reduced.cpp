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

constexpr int DIM = 3; // dimensions of the problem (2D or 3D)

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

  initMPIType(); // to be able to use MPI_VEC for Vecs
}

inline void allMPIFinalize() {
  freeMPIType();
  MPI_Finalize();
}

void gatherAndSaveAllPositions(int mpiSize, int mpiRank, int n, int steps,
                               double dt,
                               const std::vector<Vec<DIM>> &localPositions,
                               const std::vector<Vec<DIM>> &localVelocities,
                               const std::vector<double> &masses,
                               const std::string &filetag) {
  // gather all final positions and velocities to root process
  Bodies<DIM, EmptyAttributes> final;

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
        MPI_Recv(&final.position[i], 1, MPI_VEC, r, 200 + i, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(&final.velocity[i], 1, MPI_VEC, r, 400 + i, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }
    }

    utils::saveToFile("test-MPI-reduced-" + filetag + ".out", steps, dt, final,
                      false);
  }
}

void runMPIReduced(int argc, char **argv,
                   const forces::force<DIM, EmptyAttributes> &force,
                   int outputStride) {
  int mpiSize, mpiRank;
  allMPIInit(&argc, &argv, mpiSize, mpiRank);

  // if (mpiRank == 0) {
  //  utils::generateRandomToFile<DIM>("test1.in.out", 3, 1000, 0.01, 42);
  // }

  Bodies<DIM, EmptyAttributes> bodies;

  int n;
  int steps;
  double dt;

  std::vector<double> masses;
  std::vector<EmptyAttributes> attributes;

  if (mpiRank == 0) {
    // Root process reads input, and will broadcast the data.
    utils::readFromFile("test1.in.out", steps, dt, bodies);
    masses.resize(bodies.localSize());
    std::copy_n(bodies.mass.begin(), bodies.localSize(), masses.begin());
    n = bodies.globalSize();
  }

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  attributes.resize(n);

  // assert(n % mpiSize == 0); /// TODO later remove
  if (n % mpiSize != 0) {
    std::cout << "n should be divisible by mpiSize" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int locN = n / mpiSize;

  std::cout << " rank=" << mpiRank << " locN: " << locN << std::endl;

  if (mpiRank != 0) {
    masses.resize(n);
    bodies.resize(locN, locN, 0);
  }

  MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // send local positions and local velocities
  if (mpiRank == 0) {
    /* scatter cyclically: particle i goes to rank (i % comm_sz), local index =
     * i / comm_sz */
    for (int i = 0; i < n; i++) {
      int owner = i % mpiSize;
      int lidx = i / mpiSize;
      if (owner == 0) {
        bodies.position[lidx] = bodies.position[i];
        bodies.velocity[lidx] = bodies.velocity[i];
      } else {
        /* send mass and pos/vel to owner */
        MPI_Send(&bodies.position[i], 1, MPI_VEC, owner, 1, MPI_COMM_WORLD);
        MPI_Send(&bodies.velocity[i], 1, MPI_VEC, owner, 2, MPI_COMM_WORLD);
      }
    }
    bodies.resize(locN, locN, 0);
  } else {
    // get local positions and local velocities from root process
    /* receive cyclicly-scattered items from root */
    for (int i = mpiRank; i < n; i += mpiSize) {
      int lidx = i / mpiSize;
      MPI_Recv(&bodies.position[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&bodies.velocity[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }

  MPIAccumulatorReduced<DIM, EmptyAttributes> accumulator(
      MPI_VEC, n / mpiSize, mpiSize, mpiRank, masses, force, attributes);
  // integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
  // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
  integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

  for (int step = 0; step < steps; step++) {
    integrator.step(bodies, dt);

    if (step % outputStride == 0)
      gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt, bodies.position,
                                bodies.velocity, masses, std::to_string(step/outputStride));
  }

  gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt, bodies.position,
                            bodies.velocity, masses, "final");

  allMPIFinalize();
}

void run3Charges(int argc, char **argv, int outputStride) {
  int mpiSize, mpiRank;
  allMPIInit(&argc, &argv, mpiSize, mpiRank);

  Bodies<DIM, forces::charge> bodies;

  int n;
  int steps;
  double dt;

  std::vector<double> masses;
  std::vector<forces::charge> attributes;

  if (mpiRank == 0) {
    // Root process reads input, and will broadcast the data.
    steps = 1000;
    dt = 0.01;

    bodies.resize(3, 3, 0);
    bodies.local(0).pos()[0] = 1.0; bodies.local(0).pos()[1] = 1.0; bodies.local(0).pos()[2] = 1.0;
    bodies.local(0).vel()[0] = 0.0; bodies.local(0).vel()[1] = 0.0; bodies.local(0).vel()[2] = 0.0;
    bodies.local(0).mass() = 1.0;
    bodies.local(0).attributes().charge = 1e-4;

    bodies.local(1).pos()[0] = -1.0; bodies.local(1).pos()[1] = -1.0; bodies.local(1).pos()[2] = 1.0;
    bodies.local(1).vel()[0] = 0.0; bodies.local(1).vel()[1] = 0.0; bodies.local(1).vel()[2] = 0.0;
    bodies.local(1).mass() = 1.0;
    bodies.local(1).attributes().charge = 1e-4;

    bodies.local(2).pos()[0] = -1.0; bodies.local(2).pos()[1] = 1.0; bodies.local(2).pos()[2] = 1.0;
    bodies.local(2).vel()[0] = 0.0; bodies.local(2).vel()[1] = 0.0; bodies.local(2).vel()[2] = 0.0;
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

  // assert(n % mpiSize == 0); /// TODO later remove
  if (n % mpiSize != 0) {
    std::cout << "n should be divisible by mpiSize" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int locN = n / mpiSize;

  std::cout << " rank=" << mpiRank << " locN: " << locN << std::endl;

  if (mpiRank != 0) {
    masses.resize(n);
    attributes.resize(n);
    bodies.resize(locN, locN, 0);
  }

  MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(attributes.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // send local positions and local velocities
  if (mpiRank == 0) {
    /* scatter cyclically: particle i goes to rank (i % comm_sz), local index =
     * i / comm_sz */
    for (int i = 0; i < n; i++) {
      int owner = i % mpiSize;
      int lidx = i / mpiSize;
      if (owner == 0) {
        bodies.position[lidx] = bodies.position[i];
        bodies.velocity[lidx] = bodies.velocity[i];
      } else {
        /* send mass and pos/vel to owner */
        MPI_Send(&bodies.position[i], 1, MPI_VEC, owner, 1, MPI_COMM_WORLD);
        MPI_Send(&bodies.velocity[i], 1, MPI_VEC, owner, 2, MPI_COMM_WORLD);
      }
    }
    bodies.resize(locN, locN, 0);
  } else {
    // get local positions and local velocities from root process
    /* receive cyclicly-scattered items from root */
    for (int i = mpiRank; i < n; i += mpiSize) {
      int lidx = i / mpiSize;
      MPI_Recv(&bodies.position[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Recv(&bodies.velocity[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }

  forces::coulomb<DIM> coulomb{};
  MPIAccumulatorReduced<DIM, forces::charge> accumulator(
      MPI_VEC, n / mpiSize, mpiSize, mpiRank, masses, coulomb, attributes);
  integrators::Euler<DIM, forces::charge> integrator(accumulator);
  // integrators::Sympletic<DIM, forces::charge> integrator(accumulator);
  // integrators::Verlet<DIM, forces::charge> integrator(accumulator);
  // integrators::RK4<DIM, forces::charge> integrator(accumulator);

  for (int step = 0; step < steps; step++) {
    integrator.step(bodies, dt);

    if (step % outputStride == 0)
      gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt, bodies.position,
                                bodies.velocity, masses, std::to_string(step/outputStride));
  }

  gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt, bodies.position,
                            bodies.velocity, masses, "final");

  allMPIFinalize();
}

int main(int argc, char **argv) {

  run3Charges(argc, argv, 20);
  return 0;

  forces::gravity<DIM> gravity{};
  const forces::force<DIM, EmptyAttributes> &force = gravity;

  runMPIReduced(argc, argv, force, 20);
}
