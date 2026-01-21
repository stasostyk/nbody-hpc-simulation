#include "../nbody-solver.hpp"
#include "mpi-accumulator-reduced.hpp"

template <int DIM, typename Attributes>
class MPIReducedSolver : public NBodySolver<DIM, Attributes>
{
public:
    using NBodySolver<DIM, Attributes>::NBodySolver;  // Inherit constructors

    AccelerationAccumulator<DIM, Attributes> &getAccumulator() override
    {
        this->accumulator = std::make_unique<MPIAccumulatorReduced<DIM, Attributes>>(
            MPI_VEC, (this->n) / mpiSize, mpiSize, mpiRank, masses, this->force, attributes
        );
        return *(this->accumulator);
    }

    void init(int argc, char **argv, const std::string &filename) override {
        allMPIInit(&argc, &argv);
        std::cout << "MPI RANK: " << mpiRank << " (of " << mpiSize << ")" << std::endl;

        if (mpiRank == 0) {
            // Root process reads input, and will broadcast the data.
            this->readFromFile(filename);

            masses.resize((this->bodies).localSize());
            std::copy_n((this->bodies).mass.begin(), (this->bodies).localSize(), masses.begin());
        }

        MPI_Bcast(&(this->n), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(this->steps), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(this->dt), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        attributes.resize(this->n);

        if ((this->n) % mpiSize != 0) {
            std::cout << "n should be divisible by mpiSize" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        locN = (this->n) / mpiSize;

        std::cout << " rank=" << mpiRank << " locN: " << locN << std::endl;

        if (mpiRank != 0) {
            masses.resize(this->n);
            (this->bodies).resize(locN, locN, 0);
        }

        MPI_Bcast(masses.data(), this->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // send local positions and local velocities
        if (mpiRank == 0) {
            // scatter cyclically: particle i goes to rank (i % comm_sz), local index = i / comm_sz 
            for (int i = 0; i < this->n; i++) {
            int owner = i % mpiSize;
            int lidx = i / mpiSize;
            if (owner == 0) {
                (this->bodies).position[lidx] = (this->bodies).position[i];
                (this->bodies).velocity[lidx] = (this->bodies).velocity[i];
            } else {
                // send mass and pos/vel to owner
                MPI_Send(&(this->bodies).position[i], 1, MPI_VEC, owner, 1, MPI_COMM_WORLD);
                MPI_Send(&(this->bodies).velocity[i], 1, MPI_VEC, owner, 2, MPI_COMM_WORLD);
            }
            }
            (this->bodies).resize(locN, locN, 0);
        } else {
            // get local positions and local velocities from root process
            // receive cyclicly-scattered items from root 
            for (int i = mpiRank; i < (this->n); i += mpiSize) {
            int lidx = i / mpiSize;
            MPI_Recv(&(this->bodies).position[lidx], 1, MPI_VEC, 0, 1, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
            MPI_Recv(&(this->bodies).velocity[lidx], 1, MPI_VEC, 0, 2, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
            }
        }
    }

    void finalize() override {
        allMPIFinalize(); 
    }

private:
    static MPI_Datatype MPI_VEC;

    int mpiSize;
    int mpiRank;
    int locN;

    std::vector<double> masses;
    std::vector<EmptyAttributes> attributes;

    static void initMPIType() {
        MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
        MPI_Type_commit(&MPI_VEC);
    }

    static void freeMPIType() {
        MPI_Type_free(&MPI_VEC);
    }

    void allMPIFinalize() {
        freeMPIType();
        MPI_Finalize();
    }

    void allMPIInit(int *argc, char ***argv) {
        MPI_Init(argc, argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        initMPIType(); // to be able to use MPI_VEC for Vecs
    }

    void saveOutput(const std::string &filename) override {
        // gather all final positions and velocities to root process
        Bodies<DIM, EmptyAttributes> final;

        if (mpiRank == 0) {
            final.resize(this->n, this->n, 0);
            std::copy_n(masses.begin(), this->n, final.mass.begin());
        }

        if (mpiRank != 0) {
            for (int i = mpiRank, lidx = 0; i < (this->n); i += mpiSize, lidx++) {
                MPI_Send(&(this->bodies).position[lidx], 1, MPI_VEC, 0, 200 + i, MPI_COMM_WORLD);
                MPI_Send(&(this->bodies).velocity[lidx], 1, MPI_VEC, 0, 400 + i, MPI_COMM_WORLD);
            }
        }

        if (mpiRank == 0) {
            for (int i = mpiRank, lidx = 0; i < (this->n); i += mpiSize, lidx++) {
                if (mpiRank == 0) {
                    final.position[i] = (this->bodies).position[lidx];
                    final.velocity[i] = (this->bodies).velocity[lidx];
                } 
            }

            for (int r = 1; r < mpiSize; r++) {
                for (int i = r; i < this->n; i += mpiSize) {
                    MPI_Recv(&final.position[i], 1, MPI_VEC, r, 200 + i, MPI_COMM_WORLD,
                            MPI_STATUS_IGNORE);
                    MPI_Recv(&final.velocity[i], 1, MPI_VEC, r, 400 + i, MPI_COMM_WORLD,
                            MPI_STATUS_IGNORE);
                }
            }

            utils::saveToFile<DIM>(filename, this->steps, this->dt, final, false);
        }
    }

};

// Define a default value, but the actual one is initialized in mpiInit()
template <int DIM, typename Attributes>
MPI_Datatype MPIReducedSolver<DIM, Attributes>::MPI_VEC = MPI_DATATYPE_NULL;
