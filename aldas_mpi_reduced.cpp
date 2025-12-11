/*

compilation:
mpic++ -O3 -std=c++11 -o nbody_mpi_reduced_packed nbody_mpi_reduced_packed.cpp

running:
mpirun -n 4 ./nbody_mpi_reduced_packed 400 100 0.01 0

*/

#include <mpi.h>
#include <assert.h>
#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>

const static int DIM = 2; // dimensions of the problem (2D or 3D)
const double G = 6.673e-11; // gravitational constant

// using Vec = std::array<double, DIM>; // vector in our DIM-dimensional space

struct Vec : std::array<double, DIM> {
    using std::array<double, DIM>::array;  // Inherit constructors
    
    // Assignment from scalar
    Vec& operator=(double scalar) {
        std::fill(this->begin(), this->end(), scalar);
        return *this;
    }

    // Addition (member function)
    Vec operator+(const Vec& other) const {
        Vec result;
        for (size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] + other[i];
        }
        return result;
    }
    
    // Compound addition
    Vec& operator+=(const Vec& other) {
        for (size_t i = 0; i < DIM; ++i) {
            (*this)[i] += other[i];
        }
        return *this;
    }

    // Subtraction (member function)
    Vec operator-(const Vec& other) const {
        Vec result;
        for (size_t i = 0; i < DIM; ++i) {
            result[i] = (*this)[i] - other[i];
        }
        return result;
    }
    
    // Compound subtraction
    Vec& operator-=(const Vec& other) {
        for (size_t i = 0; i < DIM; ++i) {
            (*this)[i] -= other[i];
        }
        return *this;
    }

    // Compound scalar multiplication
    Vec& operator*=(double scalar) {
        for (size_t i = 0; i < DIM; ++i) {
            (*this)[i] *= scalar;
        }
        return *this;
    }

    // Scalar multiplication (member function)
    Vec operator*(double scalar) const {
        Vec result = *this;
        result *= scalar;
        return result;
    }
};

// Non-member function for scalar * Vec
Vec operator*(double scalar, const Vec& v) {
    return v * scalar;  // Reuse the member function
}

static MPI_Datatype MPI_VEC;

inline static void initMPIType() {  
    MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
    MPI_Type_commit(&MPI_VEC);
}

inline static void freeMPIType() {
    MPI_Type_free(&MPI_VEC);
}

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


inline void compareOutputsSingleFiles(std::string filename_serial, std::string filename_mpi) {
    double e = 1e-10;

    int n_serial, n_mpi;
    int steps_serial, steps_mpi;
    double dt_serial, dt_mpi;

    std::vector<double> masses_serial, masses_mpi;
    std::vector<Vec> positions_serial, positions_mpi;

    readFromFile(filename_serial, n_serial, steps_serial, dt_serial, 
                    masses_serial, positions_serial, positions_serial, false);
    readFromFile(filename_mpi, n_mpi, steps_mpi, dt_mpi, 
                    masses_mpi, positions_mpi, positions_mpi, false);

    assert(steps_mpi == steps_serial);
    assert(n_mpi == n_serial);
    assert(fabs(dt_mpi - dt_serial) < e);
    
    for (int i = 0; i < n_mpi; i++) {
        assert(fabs(masses_mpi[i]- masses_serial[i]) < e);
        for (int j = 0; j < DIM; j++) {
            assert(fabs(positions_mpi[i][j] - positions_serial[i][j]) < e);
        }
    }

    std::cout << "CHECK FINISHED, EVERYTHING IS FINE" << std::endl;

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

void runMPISerial(int argc, char** argv) {
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

// TODO assumes n%mpiSize==0, later change
// note: cyclic distribution
static inline int global_to_local(int glb, int comm_sz) {
    return glb / comm_sz;
}


/* First_index(glb_part1, owner, comm_sz, n)
 * Return the first global index assigned to 'owner' that is strictly > glb_part1.
 * (Used to find which particles in the received block to pair with.)
 *
 * In cyclic distribution owner has indices: owner, owner+P, owner+2P, ...
 * We find smallest glb_part2 = owner + k*P > glb_part1.
 */
// TODO this assumes n%mpiSize == 0
int first_index_gt(int glb_part1, int owner, int comm_sz, int n) {
    int glb = owner;
    if (glb <= glb_part1) {
        int k = (glb_part1 - owner) / comm_sz + 1;
        glb = owner + k * comm_sz;
    }
    return (glb < n) ? glb : -1;
}

void runMPIReduced(int argc, char** argv) {
    int mpiSize, mpiRank;
    allMPIInit(&argc, &argv, mpiSize, mpiRank);

    // TODO don't assume n is divisible by mpiSize

    int n;
    int steps;
    double dt;

    std::vector<double> masses;
    std::vector<Vec> allPositions;
    std::vector<Vec> allVelocities; // only filled by mpiRank=0
    std::vector<Vec> localVelocities;
    std::vector<Vec> localForces;
    std::vector<Vec> localPositions;

    if (mpiRank == 0) {
        // Root process reads input, and will broadcast the data.
        readFromFile("test1.in.out", n, steps, dt, masses, allPositions, allVelocities);
    } 

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // std::cout << " RANK r=" << mpiRank << "   is here 1" << std::endl;

    // assert(n % mpiSize == 0); /// TODO later remove
    if (n%mpiSize != 0) {
        std::cout << "n should be divisible by mpiSize" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // TODO make this computation only in root process (mpiRank 0)
    std::vector<int> counts, displs;
    make_counts_displs(n, mpiSize, counts, displs);
    int locN = counts[mpiRank];
    int locOffset = displs[mpiRank];

    std::cout << " rank=" << mpiRank << " locN: " << locN << " locOffset: " << locOffset << std::endl;

    if (mpiRank != 0) {
        masses.resize(n);
        // positions.resize(n);
    }

    localPositions.resize(locN);
    localVelocities.resize(locN);
    localForces.resize(locN);

    MPI_Bcast(masses.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // send local positions and local velocities
    if (mpiRank == 0) {
        /* scatter cyclically: particle i goes to rank (i % comm_sz), local index = i / comm_sz */
        for (int i = 0; i < n; i++) {
            int owner = i % mpiSize;
            int lidx = i / mpiSize;
            if (owner == 0) {
                localPositions[lidx] = allPositions[i];
                localVelocities[lidx] = allVelocities[i];
            } else {
                /* send mass and pos/vel to owner */
                MPI_Send(&allPositions[i], 1, MPI_VEC, owner, 1, MPI_COMM_WORLD);
                MPI_Send(&allVelocities[i], 1, MPI_VEC, owner, 2, MPI_COMM_WORLD);
            }
        }
    } else {
        // get local positions and local velocities from root process
        /* receive cyclicly-scattered items from root */
        for (int i = mpiRank; i < n; i += mpiSize) {
            int lidx = i / mpiSize;
            MPI_Recv(&localPositions[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&localVelocities[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    Vec *tmpData = (Vec *) malloc(sizeof(Vec) * 2 * locN);
    Vec *tmpPos = tmpData;
    Vec *tmpForces = tmpData + locN;

    int source = (mpiRank + 1) % mpiSize;
    int dest = (mpiRank - 1 + mpiSize) % mpiSize;

    // MPI_Barrier(MPI_COMM_WORLD); // TODO is this needed???

    for (int step = 0; step < steps; step++) {
        // make forces equal to zero and prepare tmpData
        for (int i = 0; i < locN; i++) {
            localForces[i] = 0.;
            tmpPos[i] = localPositions[i];
            tmpForces[i] = 0.;
        }

        // first, compute forces due to interactions among local particles
        for (int locI = 0; locI < locN; locI++) {
            int glbI = mpiRank + locI * mpiSize;
            for (int locJ = locI + 1; locJ < locN; locJ++) {
                int glbJ = mpiRank + locJ * mpiSize;
                Vec force = force_on_p1(localPositions[locI], localPositions[locJ],
                                        masses[glbI], masses[glbJ]);
                
                localForces[locI] += force;
                tmpForces[locJ] -= force;
            }
        }

        // ring pass: mpiSize-1 phases
        for (int phase = 1; phase < mpiSize; phase++) {
            // send tmpData (tmp Positions and tmp Forces) to dest, and receive into the same 
            // buffer from the source
            // TODO why the tag is 99??
            MPI_Sendrecv_replace(
                tmpData, 2*locN, MPI_VEC,
                dest, 99,
                source, 99,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );

            int owner = (mpiRank + phase) % mpiSize;

            // compute interactions between local particles and the received ones
            for (int locI = 0, glbI = mpiRank; locI < locN; locI++, glbI += mpiSize) {
                // find first global index > glbI that belongs to owner
                int glbJ = first_index_gt(glbI, owner, mpiSize, n);
                
                while (glbJ != -1 && glbJ < glbI + mpiSize && glbJ < n) {
                    int locJ = global_to_local(glbJ, mpiSize);
                    
                    // glbI < glbJ by construction
                    Vec force = force_on_p1(
                        localPositions[locI], localPositions[locJ],
                        masses[glbI], masses[glbJ]
                    );

                    localForces[locI] += force;
                    tmpForces[locJ] -= force;

                    glbJ += mpiSize;
                    if (glbJ >= n) break;
                }
            }
        }

        // one last send/recv of tmp arrays
        // TODO why tag is 98??
        MPI_Sendrecv_replace(
            tmpData, 2*locN, MPI_VEC,
            dest, 98, source, 98,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );

        // add the tmpForces that were just received
        for (int i = 0; i < locN; i++) {
            localForces[i] += tmpForces[i];
        }

        // update velocities and positions (Euler)
        for (int locI = 0, glbI = mpiRank; locI < locN; locI++, glbI += mpiSize) {
            double invM = 1. / masses[glbI];
            localPositions[locI] += dt * localVelocities[locI];
            localVelocities[locI] += dt * invM * localForces[locI];
        }
    }

    // gather all final positions and velocities to root process
    std::vector<Vec> finalPositions;
    std::vector<Vec> finalVelocities;

    if (mpiRank == 0) {
        finalPositions.resize(n);
        finalVelocities.resize(n);
    }

    for (int i = mpiRank, lidx = 0; i < n; i += mpiSize, lidx++) {
        if (mpiRank == 0) {
            finalPositions[i] = localPositions[lidx];
            finalVelocities[i] = localVelocities[lidx];
        } else {
            MPI_Send(&localPositions[lidx], 1, MPI_VEC, 0, 200+i, MPI_COMM_WORLD);
            MPI_Send(&localVelocities[lidx], 1, MPI_VEC, 0, 400+i, MPI_COMM_WORLD);
        }
    }

    if (mpiRank == 0) {
        for (int r = 1; r < mpiSize; r++) {
            for (int i = r; i < n; i += mpiSize) {
                MPI_Recv(&finalPositions[i], 1, MPI_VEC, r, 200+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&finalVelocities[i], 1, MPI_VEC, r, 400+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        saveToFile("test-MPI-reduced-final.out", n, steps, dt, masses, finalPositions, finalVelocities, false);
    }

    free(tmpData);
    allMPIFinalize();
}

int main(int argc, char** argv) {

    // runSerial();
    // return 0;

    // compareOutputs();
    // return 0;

    // runMPISerial();

    runMPIReduced(argc, argv);
    // compareOutputsSingleFiles("test1.9.out", "test-MPI-reduced-final.out");
}