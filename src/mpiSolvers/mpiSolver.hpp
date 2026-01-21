#include "../nbody-solver.hpp"
#include "mpi-accumulator.hpp"

template <int DIM, typename Attributes>
class MPISolver : public NBodySolver<DIM, Attributes>
{
public:
    using NBodySolver<DIM, Attributes>::NBodySolver;  // Inherit constructors

    AccelerationAccumulator<DIM, Attributes> &getAccumulator() override
    {
        this->accumulator = std::make_unique<MPIAccumulator<DIM, Attributes>>(
            MPI_VEC, locN, counts, displs, this->force
        );
        return *(this->accumulator);
    }

    void init(int argc, char **argv, const std::string &filename) override {
        allMPIInit(&argc, &argv);
        std::cout << "MPI RANK: " << mpiRank << " (of " << mpiSize << ")" << std::endl;
        if (mpiRank == 0) {
            this->readFromFile(filename);
        }

        MPI_Bcast(&(this->n), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(this->steps), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(this->dt), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        make_counts_displs(this->n, mpiSize);
        locN = counts[mpiRank];
        locOffset = displs[mpiRank];

        (this->bodies).resize(this->n, locN, locOffset);

        MPI_Bcast((this->bodies).position.data(), this->n, MPI_VEC, 0, MPI_COMM_WORLD);
        MPI_Bcast((this->bodies).velocity.data(), this->n, MPI_VEC, 0, MPI_COMM_WORLD);
        MPI_Bcast((this->bodies).mass.data(), this->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void finalize() override {
        allMPIFinalize();
    }

private:
    static MPI_Datatype MPI_VEC;

    int mpiSize;
    int mpiRank;
    int locN;
    int locOffset;
    std::vector<int> counts;
    std::vector<int> displs;

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

    bool canPrint() override {
        return mpiRank == 0;
    }

    // utility: block counts/displacements for an array of N elements across comm_sz
    void make_counts_displs(int N, int comm_sz) {
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

    void oneStep() override { // after integrator.step()
        MPI_Allgatherv((this->bodies).position.data() + (this->bodies).localOffset(),
                        (this->bodies).localSize(), MPI_VEC, (this->bodies).position.data(),
                        counts.data(), displs.data(), MPI_VEC, MPI_COMM_WORLD);
    }

};


// Define a default value, but the actual one is initialized in mpiInit()
template <int DIM, typename Attributes>
MPI_Datatype MPISolver<DIM, Attributes>::MPI_VEC = MPI_DATATYPE_NULL;
