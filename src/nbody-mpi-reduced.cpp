#include <mpi.h>
#include <assert.h>
#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "utils.hpp"
#include "forces/gravity.hpp"

constexpr int DIM = 2; // dimensions of the problem (2D or 3D)

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

inline void compareOutputsSingleFiles(std::string filename_serial, std::string filename_mpi) {
    double e = 1e-2;

    int n_serial, n_mpi;
    int steps_serial, steps_mpi;
    double dt_serial, dt_mpi;

    std::vector<double> masses_serial, masses_mpi;
    std::vector<Vec<DIM>> positions_serial, positions_mpi;

    utils::readFromFile(filename_serial, n_serial, steps_serial, dt_serial, 
                    masses_serial, positions_serial, positions_serial, false);
    utils::readFromFile(filename_mpi, n_mpi, steps_mpi, dt_mpi, 
                    masses_mpi, positions_mpi, positions_mpi, false);

    assert(steps_mpi == steps_serial);
    assert(n_mpi == n_serial);
    assert(fabs(dt_mpi - dt_serial) < e);
    
    for (int i = 0; i < n_mpi; i++) {
        assert(fabs(masses_mpi[i]- masses_serial[i]) < e);
        for (int j = 0; j < DIM; j++) {
            assert(fabs(positions_mpi[i][j] - positions_serial[i][j]) < e);
            if (fabs(positions_mpi[i][j] - positions_serial[i][j]) > e) {
                std::cout << "PARTICLE i=" << i << " IS WRONG" << std::endl;
                std::cout << "   diff: " << positions_mpi[i][j] << " and " << positions_serial[i][j] << std::endl;
            }
        }
    }

    std::cout << "CHECK FINISHED, EVERYTHING IS FINE" << std::endl;

}


// assumes n%mpiSize==0, cyclic distribution
static inline int global_to_local(int glb, int comm_sz) {
    return glb / comm_sz;
}

/* First_index(glb_part1, owner, comm_sz, n)
 * Return the first global index assigned to 'owner' that is strictly > glb_part1.
 * (Used to find which particles in the received block to pair with.)
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

void runMPIReduced(int argc, char** argv, const forces::force<DIM>& force) {
    int mpiSize, mpiRank;
    allMPIInit(&argc, &argv, mpiSize, mpiRank);

    int n;
    int steps;
    double dt;

    std::vector<double> masses;
    std::vector<Vec<DIM>> allPositions;
    std::vector<Vec<DIM>> allVelocities; // only filled by mpiRank=0
    std::vector<Vec<DIM>> localVelocities;
    std::vector<Vec<DIM>> localForces;
    std::vector<Vec<DIM>> localPositions;

    if (mpiRank == 0) {
        // Root process reads input, and will broadcast the data.
        utils::readFromFile("test1.in.out", n, steps, dt, masses, allPositions, allVelocities);
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

    int locN = n / mpiSize;

    std::cout << " rank=" << mpiRank << " locN: " << locN << std::endl;

    if (mpiRank != 0) {
        masses.resize(n);
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

    Vec<DIM> *tmpData = (Vec<DIM> *) malloc(sizeof(Vec<DIM>) * 2 * locN);
    Vec<DIM> *tmpPos = tmpData;
    Vec<DIM> *tmpForces = tmpData + locN;

    int source = (mpiRank + 1) % mpiSize;
    int dest = (mpiRank - 1 + mpiSize) % mpiSize;

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
                Vec<DIM> f = force(localPositions[locI], masses[glbI], 
                                   localPositions[locJ], masses[glbJ]);
                
                localForces[locI] += f;
                tmpForces[locJ] -= f;
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
                    Vec<DIM> f = force(
                        localPositions[locI], masses[glbI],
                        tmpPos[locJ], masses[glbJ]
                    );

                    localForces[locI] += f;
                    tmpForces[locJ] -= f;

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
    std::vector<Vec<DIM>> finalPositions;
    std::vector<Vec<DIM>> finalVelocities;

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

        utils::saveToFile("test-MPI-reduced-final.out", n, steps, dt, masses, finalPositions, finalVelocities, false);
    }

    free(tmpData);
    allMPIFinalize();
}

int main(int argc, char** argv) {
    
    forces::gravity<DIM> gravity{};
    const forces::force<DIM> &force = gravity;

    runMPIReduced(argc, argv, force);
    // compareOutputsSingleFiles("../build/test1-MPI.999.out", "../build/test-MPI-reduced-final.out");
}