#include <algorithm>
#include <mpi.h>
#include <assert.h>
#include <vector>
#include <cmath>
#include <iostream>

#include "integrators/euler.hpp"
#include "integrators/sympletic.hpp"
#include "integrators/verlet.hpp"
#include "integrators/rk4.hpp"
#include "utils.hpp"
#include "forces/gravity.hpp"
#include "mpi-accumulator-reduced.hpp"

constexpr int DIM = 3; // dimensions of the problem (2D or 3D)

static MPI_Datatype MPI_VEC;

inline static void initMPIType() {  
    MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
    MPI_Type_commit(&MPI_VEC);
}

inline static void freeMPIType() {
    MPI_Type_free(&MPI_VEC);
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

void gatherAndSaveAllPositions(int mpiSize, int mpiRank, int n, int steps, double dt,
                                const std::vector<Vec<DIM>> &localPositions,
                                const std::vector<Vec<DIM>> &localVelocities,
                                const std::vector<double> &masses,
                                const std::string &filetag) {
    // gather all final positions and velocities to root process
    bodies<DIM> final;

    if (mpiRank == 0) {
        final.resize(n, n, 0);
    }

    for (int i = mpiRank, lidx = 0; i < n; i += mpiSize, lidx++) {
        if (mpiRank == 0) {
            final.position[i] = localPositions[lidx];
            final.velocity[i] = localVelocities[lidx];
        } else {
            MPI_Send(&localPositions[lidx], 1, MPI_VEC, 0, 200+i, MPI_COMM_WORLD);
            MPI_Send(&localVelocities[lidx], 1, MPI_VEC, 0, 400+i, MPI_COMM_WORLD);
        }
    }

    if (mpiRank == 0) {
        std::copy_n(masses.begin(), n, final.mass.begin());
        for (int r = 1; r < mpiSize; r++) {
            for (int i = r; i < n; i += mpiSize) {
                MPI_Recv(&final.position[i], 1, MPI_VEC, r, 200+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&final.velocity[i], 1, MPI_VEC, r, 400+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        utils::saveToFile("test-MPI-reduced-" + filetag + ".out", steps, dt, final, false);
    }
}

void runMPIReduced(int argc, char** argv, const forces::force<DIM>& force, int outputStride) {
    int mpiSize, mpiRank;
    allMPIInit(&argc, &argv, mpiSize, mpiRank);

    if (mpiRank == 0) {
        utils::generateRandomToFile<DIM>("test1.in.out", 3, 1000, 0.01, 42);
    }

    bodies<DIM> bodies;

    int n;
    int steps;
    double dt;

    std::vector<double> masses;

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

    // assert(n % mpiSize == 0); /// TODO later remove
    if (n%mpiSize != 0) {
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
        /* scatter cyclically: particle i goes to rank (i % comm_sz), local index = i / comm_sz */
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
            MPI_Recv(&bodies.position[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&bodies.velocity[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    MPIAccumulatorReduced<DIM> accumulator(MPI_VEC, n, mpiSize, mpiRank, masses, force);
    //integrators::Euler<DIM> integrator(accumulator);
    //integrators::Sympletic<DIM> integrator(accumulator);
    //integrators::Verlet<DIM> integrator(accumulator);
    integrators::RK4<DIM> integrator(accumulator);

    for (int step = 0; step < steps; step++) {
        integrator.step(bodies, dt);

        if (step % outputStride == 0)
            gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt, bodies.position, bodies.velocity, masses, std::to_string(step));
    }

    gatherAndSaveAllPositions(mpiSize, mpiRank, n, steps, dt, bodies.position, bodies.velocity, masses, "final");

    allMPIFinalize();
}

int main(int argc, char** argv) {
    
    forces::gravity<DIM> gravity{};
    const forces::force<DIM> &force = gravity;

    runMPIReduced(argc, argv, force, 20);
}
