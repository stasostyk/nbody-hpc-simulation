#include <mpi.h>
#include <assert.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "forces/gravity.hpp"
#include "utils.hpp"

constexpr int DIM = 2;
constexpr double EPS = 1e-3;

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

double critical_timestep(double mass, const Vec<DIM> &vel, const Vec<DIM> &force, double dt_max) {
    double acc_norm_sq = 0.0;
    double vel_norm_sq = 0.0;
    for (int d = 0; d < DIM; ++d) {
        double acc = force[d] / mass;
        acc_norm_sq += acc * acc;
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

double update_timesteps(const std::vector<Vec<DIM>> &local_vel,
                        const std::vector<Vec<DIM>> &local_forces,
                        const std::vector<double> &masses,
                        int locOffset,
                        double dt_max,
                        std::vector<double> &particle_dt) {
    const double min_dt = dt_max * 1e-6;
    double next_dt = dt_max;

    for (size_t i = 0; i < local_vel.size(); ++i) {
        int q = locOffset + static_cast<int>(i);
        double dt_crit = critical_timestep(masses[q], local_vel[i], local_forces[i], dt_max);
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

inline void runSerial(const forces::force<DIM>& force) {

    // BELOW CODE: for generating random test1.in.out
    // utils::generateRandomToFile<DIM>("test1.in.out");
    // exit(0);

    int n;
    double total_time;
    double dt_max;

    std::vector<double> masses;
    std::vector<Vec<DIM>> positions;
    std::vector<Vec<DIM>> velocities;
    std::vector<Vec<DIM>> forces;
    std::vector<double> particle_dt;

    utils::readFromFile("test1.in.out", n, total_time, dt_max, masses, positions, velocities);

    forces.resize(n);
    particle_dt.assign(n, dt_max);
    for (int i = 0; i < n; i++) for (int j = 0; j < DIM; j++)
        forces[i][j] = 0.;

    double time = 0.0;
    int step = 0;
    while (time < total_time) {
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
                forces[i] += force(positions[i], masses[i], positions[k], masses[k]);
            }
        }

        double next_dt = update_timesteps(velocities, forces, masses, 0, dt_max, particle_dt);
        double dt = std::min(next_dt, total_time - time);
        if (dt <= 0.0) {
            break;
        }

        // apply forces, i.e. calc new positions and velocities
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < DIM; j++) {
                positions[i][j] += dt * velocities[i][j];
                velocities[i][j] += dt / masses[i] * forces[i][j];
            }
        }
    
        // save for testing
        utils::saveToFile("test1." + std::to_string(step) + ".out", n, total_time, dt_max,
            masses, positions, velocities, false);
        time += dt;
        ++step;
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
        double total_time_serial, total_time_mpi;
        double dt_max_serial, dt_max_mpi;

        std::vector<double> masses_serial, masses_mpi;
        std::vector<Vec<DIM>> positions_serial, positions_mpi;

        std::string filename_suf = + "." + std::to_string(step) + ".out";
        std::string filename_serial = serial_pref + filename_suf;
        std::string filename_mpi = mpi_pref + filename_suf;

        utils::readFromFile(filename_serial, n_serial, total_time_serial, dt_max_serial, 
                     masses_serial, positions_serial, positions_serial, false);
        utils::readFromFile(filename_mpi, n_mpi, total_time_mpi, dt_max_mpi, 
                     masses_mpi, positions_mpi, positions_mpi, false);

        assert(n_mpi == n_serial);
        assert(fabs(total_time_mpi - total_time_serial) < e);
        assert(fabs(dt_max_mpi - dt_max_serial) < e);
        
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

    forces::gravity<DIM> gravity{};
    const forces::force<DIM> &force = gravity;

    int mpiSize, mpiRank;
    allMPIInit(&argc, &argv, mpiSize, mpiRank);
    
    int n;
    double total_time;
    double dt_max;

    std::vector<double> masses;
    std::vector<Vec<DIM>> positions;
    std::vector<Vec<DIM>> allVelocities; // only filled by mpiRank=0
    std::vector<Vec<DIM>> localVelocities;
    std::vector<Vec<DIM>> localForces;

    if (mpiRank == 0) {
        // Root process reads input, and will broadcast the data.
        utils::readFromFile("test1.in.out", n, total_time, dt_max, masses, positions, allVelocities);
    } 

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // std::cout << " RANK r=" << mpiRank << "   is here 1" << std::endl;

    // TODO make this computation only in root process (mpiRank 0)
    std::vector<int> counts, displs;
    make_counts_displs(n, mpiSize, counts, displs);
    int locN = counts[mpiRank];
    int locOffset = displs[mpiRank];

    std::cout << " rank=" << mpiRank << " locN: " << locN << " locOffset: " << locOffset << std::endl;

    if (mpiRank != 0) {
        masses.resize(n);
        positions.resize(n);
    }

    localVelocities.resize(locN);
    localForces.resize(locN);
    std::vector<double> particle_dt(locN, dt_max);

    MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(positions.data(), n, MPI_VEC, 0, MPI_COMM_WORLD);

    // TODO broadcast velocities in place ???
    std::vector<int> sendcounts(mpiSize), senddispls(mpiSize);
    for (int p = 0; p < mpiSize; ++p) { sendcounts[p] = counts[p]; senddispls[p] = displs[p]; }
    MPI_Scatterv(allVelocities.data(), sendcounts.data(), senddispls.data(), MPI_VEC, 
                 localVelocities.data(), locN, MPI_VEC, 0, MPI_COMM_WORLD);

    // std::cout << " RANK r=" << mpiRank << "   is here 2" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD); // TODO is this needed???

    double time = 0.0;
    int step = 0;

    while (time < total_time) {
        // make forces equal to zero
        for (int i = 0; i < locN; i++) {
            for (int j = 0; j < DIM; j++) {
                localForces[i][j] = 0.;
            }
        }

        // compute all forces
        for (int i = 0; i < locN; i++) {
            int q = locOffset + i;
            for (int k = 0; k < n; k++) {
                if (q == k) continue;
                localForces[i] += force(positions[q], masses[q], positions[k], masses[k]);
            }
        }

        double local_next_dt = update_timesteps(localVelocities, localForces, masses, locOffset, dt_max, particle_dt);
        double next_dt = dt_max;
        MPI_Allreduce(&local_next_dt, &next_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        double dt = std::min(next_dt, total_time - time);
        if (dt <= 0.0) {
            break;
        }

        // apply forces, i.e. calc new positions and velocities
        for (int i = 0; i < locN; i++) {
            int q = locOffset + i;
            for (int j = 0; j < DIM; j++) {
                // TODO first update velocities, and then positions, or vice versa?
                positions[q][j] += dt * localVelocities[i][j];
                localVelocities[i][j] += dt / masses[q] * localForces[i][j];
            }
        }
    
        // All tasks gather all positions (will be required for the next step)
        std::vector<int> recvcounts(mpiSize), recvdispls(mpiSize);
        for (int p = 0; p < mpiSize; ++p) { recvcounts[p] = counts[p]; recvdispls[p] = displs[p]; }
        MPI_Allgatherv(positions.data() + locOffset, locN, MPI_VEC,
                       positions.data(), recvcounts.data(), recvdispls.data(), 
                       MPI_VEC, MPI_COMM_WORLD); 

        if (mpiRank == 0) {
            // save for testing
            utils::saveToFile("test1-MPI." + std::to_string(step) + ".out", n, total_time, dt_max,
                masses, positions, allVelocities, false);
        }

        time += dt;
        ++step;
    }

    allMPIFinalize();
}
