#include <mpi.h>
#include <assert.h>
#include <vector>
#include <cmath>
#include <iostream>

#include "acceleration-accumulator.hpp"
#include "body.hpp"
#include "forces/gravity.hpp"
#include "integrators/euler.hpp"
#include "integrators/integrator.hpp"
#include "integrators/sympletic.hpp"
#include "integrators/verlet.hpp"
#include "integrators/rk4.hpp"
#include "mpi-accumulator.hpp"
#include "utils.hpp"

constexpr int DIM = 2;

static MPI_Datatype MPI_VEC;

inline static void initMPIType() {
    MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
    MPI_Type_commit(&MPI_VEC);
}

inline static void freeMPIType() {
    MPI_Type_free(&MPI_VEC);
}

// struct Particle {
//     double mass;
//     double velocity;
//     std::array<double, DIM> position;
//     std::array<double, DIM> force;

//     Particle() : mass(0.), velocity(0.), position{}, force{} {}
// };

// // gets force acting on particle p1
// inline std::array<double, DIM> force_on_p1(const Particle &p1, const Particle &p2)
// {
//     std::array<double, DIM> dist;
//     double distSquared = 1e-12; // softening to avoid division by zero
//     for (int i = 0; i < DIM; i++) {
//         dist[i] = p1.position[i] - p2.position[i];
//         distSquared += dist[i] * dist[i];
//     }

//     double dist3 = distSquared * std::sqrt(distSquared);
//     double coeff = -G * p1.mass * p2.mass / dist3;

//     std::array<double, DIM> force;
//     for (int i = 0; i < DIM; i++) 
//         force[i] = dist[i] * coeff;
    
//     return force;
// }

// utility: block counts/displacements for an array of N elements across comm_sz
void make_counts_displs(int N, int comm_sz, std::vector<int> &counts, std::vector<int> &displs) {
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

inline void runSerial(const forces::force<DIM>& force) {

    // BELOW CODE: for generating random test1.in.out
    // utils::generateRandomToFile<DIM>("test1.in.out");
    // exit(0);

    int n;
    int steps;
    double dt;

    std::vector<double> masses;
    std::vector<Vec<DIM>> positions;
    std::vector<Vec<DIM>> velocities;
    std::vector<Vec<DIM>> forces;

    // utils::readFromFile("test1.in.out", n, steps, dt, masses, positions, velocities);

    forces.resize(n);
    for (int i = 0; i < n; i++) for (int j = 0; j < DIM; j++)
        forces[i][j] = 0.;

    for (int step = 0; step < steps; step++) {
        // make forces equal to zero
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < DIM; j++) {
                forces[i][j] = 0.;
            }
        }

        // compute all forces
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                if (i == k) continue;
                //forces[i] += force(positions[i], masses[i], positions[k], masses[k]);
            }
        }

        // apply forces, i.e. calc new positions and velocities
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < DIM; j++) {
                positions[i][j] += dt * velocities[i][j];
                velocities[i][j] += dt / masses[i] * forces[i][j];
            }
        }
    
        // save for testing
        // utils::saveToFile("test1." + std::to_string(step) + ".out", n, steps, dt,
        //    masses, positions, velocities, false);
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

inline void compareOutputs() {

    std::string serial_pref = "test1";
    std::string mpi_pref = "test1-MPI";
    int steps = 10;
    double e = 1e-10;

    for (int step = 0; step < steps; step++) {
        int n_serial, n_mpi;
        int steps_serial, steps_mpi;
        double dt_serial, dt_mpi;

        std::vector<double> masses_serial, masses_mpi;
        std::vector<Vec<DIM>> positions_serial, positions_mpi;

        std::string filename_suf = + "." + std::to_string(step) + ".out";
        std::string filename_serial = serial_pref + filename_suf;
        std::string filename_mpi = mpi_pref + filename_suf;

        // utils::readFromFile(filename_serial, n_serial, steps_serial, dt_serial, 
        //              masses_serial, positions_serial, positions_serial, false);
        // utils::readFromFile(filename_mpi, n_mpi, steps_mpi, dt_mpi, 
        //             masses_mpi, positions_mpi, positions_mpi, false);

        assert(steps == steps_mpi && steps_mpi == steps_serial);
        assert(n_mpi == n_serial);
        assert(fabs(dt_mpi - dt_serial) < e);
        
        for (int i = 0; i < n_mpi; i++) {
            assert(fabs(masses_mpi[i]- masses_serial[i]) < e);
            for (int j = 0; j < DIM; j++) {
                assert(fabs(positions_mpi[i][j] - positions_serial[i][j]) < e);
            }
        }
    }

    std::cout << "CHECK FINISHED, EVERYTHING IS FINE" << std::endl;

}

int main(int argc, char** argv) {

    // runSerial(forces::gravity<DIM>{});
    // return 0;

    // compareOutputs();
    // return 0;

    /*
    TODO
    - compare results with serial version
    - later: instead of basic version, do the reduced version
    */

    constexpr int outputStride = 1;

    const forces::force<DIM> &force = forces::gravity<DIM>();
    bodies<DIM> bodies;

    int mpiSize, mpiRank;
    allMPIInit(&argc, &argv, mpiSize, mpiRank);
    
    int n;
    int steps;
    double dt;

    if (mpiRank == 0) {
        // Root process reads input, and will broadcast the data.
        utils::readFromFile("test1.in.out", steps, dt, bodies);
        n = bodies.globalSize();
    } 

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
    
    MPIAccumulator<DIM> accumulator(MPI_VEC, locN, counts, displs, force);
    // integrators::Euler<DIM> integrator(accumulator);
    // integrators::Sympletic<DIM> integrator(accumulator);
    // integrators::Verlet<DIM> integrator(accumulator);
    integrators::RK4<DIM> integrator(accumulator);

    for (int step = 0; step < steps; step++) {
        integrator.step(bodies, dt);

        if (step % outputStride == 0) {
            MPI_Allgatherv(bodies.position.data() + bodies.localOffset(),
                           bodies.localSize(), MPI_VEC, bodies.position.data(),
                           counts.data(), displs.data(), MPI_VEC, MPI_COMM_WORLD);

            if (mpiRank == 0) {
                // save for testing
                utils::saveToFile("test1-MPI." + std::to_string(step) + ".out", 
                    steps, dt, bodies, false);
            }
        }
    }

    allMPIFinalize();
}
