#include <mpi.h>
#include <assert.h>
#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>

const static int DIM = 2; // dimensions of the problem (2D or 3D)
const double G = 6.673e-11; // gravitational constant

using Vec = std::array<double, DIM>; // vector in our DIM-dimensional space

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


// gets force acting on particle p1
inline Vec force_on_p1(
    const Vec &pos1, 
    const Vec &pos2,
    double m1, double m2
) {
    Vec dist;
    double distSquared = 1e-12; // softening to avoid division by zero
    for (int i = 0; i < DIM; i++) {
        dist[i] = pos1[i] - pos2[i];
        distSquared += dist[i] * dist[i];
    }

    double dist3 = distSquared * std::sqrt(distSquared);
    double coeff = -G * m1 * m2 / dist3;

    Vec force;
    for (int i = 0; i < DIM; i++) 
        force[i] = dist[i] * coeff;
    
    return force;
}

inline void saveToFile(
    const std::string &fileName,
    int n, int steps, double dt,
    const std::vector<double> &m, 
    const std::vector<Vec> &p, 
    const std::vector<Vec> &v,
    bool saveVelocities = true
) {
    std::ofstream fout(fileName);

    fout << DIM << "\n";
    fout << n << " " << steps << " " << dt << "\n";
    for (int i = 0; i < n; i++) {
        fout << m[i] << " ";
        for (int j = 0; j < DIM; j++) fout << p[i][j] << " ";
        if (saveVelocities) {
            for (int j = 0; j < DIM; j++) fout << v[i][j] << " ";
        }
        fout << "\n";
    }
    fout.flush();
    fout.close();

    std::cout << "Saved to file " << fileName << std::endl;
}

inline void readFromFile(
    const std::string &fileName,
    int &n, int &steps, double &dt,
    std::vector<double> &m, 
    std::vector<Vec> &p, 
    std::vector<Vec> &v,
    bool readVelocities = true
) {
    std::ifstream fin(fileName);

    int dimInFile;
    fin >> dimInFile;
    assert(dimInFile == DIM);

    fin >> n >> steps >> dt;

    m.resize(n);
    p.resize(n);
    v.resize(n);

    for (int i = 0; i < n; i++) {
        fin >> m[i];
        for (int j = 0; j < DIM; j++) fin >> p[i][j];
        if (readVelocities) {
            for (int j = 0; j < DIM; j++) fin >> v[i][j];
        }
    }

    fin.close();
    std::cout << "Read from file " << fileName << std::endl;
    std::cout << "  DIM=" << DIM << " n=" << n << " steps=" << steps << " dt=" << dt << std::endl;
}

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

inline void initAndSaveRandom() {
    int n = 100;
    int steps = 10;
    double dt = 0.01;

    std::vector<double> masses(n);
    std::vector<Vec> positions(n);
    std::vector<Vec> velocities(n);

    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> posd(-100, 100);
    std::uniform_real_distribution<double> massd(1, 100);
    std::uniform_real_distribution<double> veld(-10, 10);
    for (int i = 0; i < n; ++i) {
        masses[i] = massd(rng);
        for (int j = 0; j < DIM; j++) {
            positions[i][j] = posd(rng);
            velocities[i][j] = veld(rng);
        }
    }

    // save to file
    saveToFile("test1.in.out", n, steps, dt, masses, positions, velocities);
}

inline void runSerial() {

    // BELOW CODE: for generating random test1.in.out
    // initAndSaveRandom();
    // // return 0;
    // exit(0);

    int n;
    int steps;
    double dt;

    std::vector<double> masses;
    std::vector<Vec> positions;
    std::vector<Vec> velocities;
    std::vector<Vec> forces;

    readFromFile("test1.in.out", n, steps, dt, masses, positions, velocities);

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
                Vec force = force_on_p1(positions[i], positions[k], masses[i], masses[k]);
                for (int j = 0; j < DIM; j++) {
                    forces[i][j] += force[j];
                }
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
        saveToFile("test1." + std::to_string(step) + ".out", n, steps, dt,
            masses, positions, velocities, false);
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
        std::vector<Vec> positions_serial, positions_mpi;

        std::string filename_suf = + "." + std::to_string(step) + ".out";
        std::string filename_serial = serial_pref + filename_suf;
        std::string filename_mpi = mpi_pref + filename_suf;

        readFromFile(filename_serial, n_serial, steps_serial, dt_serial, 
                     masses_serial, positions_serial, positions_serial, false);
        readFromFile(filename_mpi, n_mpi, steps_mpi, dt_mpi, 
                     masses_mpi, positions_mpi, positions_mpi, false);

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

    // runSerial();
    // return 0;

    // compareOutputs();
    // return 0;

    /*
    TODO
    - compare results with serial version
    - later: instead of basic version, do the reduced version
    */

    int mpiSize, mpiRank;
    allMPIInit(&argc, &argv, mpiSize, mpiRank);
    
    int n;
    int steps;
    double dt;

    std::vector<double> masses;
    std::vector<Vec> positions;
    std::vector<Vec> allVelocities; // only filled by mpiRank=0
    std::vector<Vec> localVelocities;
    std::vector<Vec> localForces;

    if (mpiRank == 0) {
        // Root process reads input, and will broadcast the data.
        readFromFile("test1.in.out", n, steps, dt, masses, positions, allVelocities);
    } 

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

    MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(positions.data(), n, MPI_VEC, 0, MPI_COMM_WORLD);

    // TODO broadcast velocities in place ???
    std::vector<int> sendcounts(mpiSize), senddispls(mpiSize);
    for (int p = 0; p < mpiSize; ++p) { sendcounts[p] = counts[p]; senddispls[p] = displs[p]; }
    MPI_Scatterv(allVelocities.data(), sendcounts.data(), senddispls.data(), MPI_VEC, 
                 localVelocities.data(), locN, MPI_VEC, 0, MPI_COMM_WORLD);

    // std::cout << " RANK r=" << mpiRank << "   is here 2" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD); // TODO is this needed???

    for (int step = 0; step < steps; step++) {
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

                Vec force = force_on_p1(positions[q], positions[k], masses[q], masses[k]);
                for (int j = 0; j < DIM; j++) {
                    localForces[i][j] += force[j];
                }
            }
        }

        // apply forces, i.e. calc new positions and velocities
        for (int i = 0; i < locN; i++) {
            int q = locOffset + i;
            for (int j = 0; j < DIM; j++) {
                // TODO first update velocities, and then positions, or vice versa?
                positions[q][j] += dt * localVelocities[i][j];
                localVelocities[i][j] += dt / masses[i] * localForces[i][j];
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
            saveToFile("test1-MPI." + std::to_string(step) + ".out", n, steps, dt,
                masses, positions, allVelocities, false);
        }
    }

    allMPIFinalize();
}
