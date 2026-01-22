#include <memory>

#include "nbody-solver.hpp"
#include "serialSolvers/serialSolver.hpp"
#include "serialSolvers/serialReducedSolver.hpp"
#include "mpiSolvers/mpiSolver.hpp"
#include "mpiSolvers/mpiReducedSolver.hpp"
#include "barnes_hut/barnesHutSolver.hpp"
#include "acceleration-accumulator.hpp"
#include "body.hpp"
#include "utils.hpp"
#include "forces/func.hpp"
#include "forces/gravity.hpp"
#include "integrators/integrator.hpp"
#include "integrators/euler.hpp"
#include "integrators/rk4.hpp"
#include "integrators/sympletic.hpp"
#include "integrators/verlet.hpp"

#if USE_OPENMP
    #include <omp.h>
#endif

#define DIM 2

void printHelp(char **argv) {
    std::cerr << "Usage: " << argv[0]
                << " [inputFilename] [solver] [extra]\n";
    std::cerr << "  solver: Serial SerialReduced MPI MPIReduced BarnesHut\n";
    std::cerr << "  in case BarnesHut is used, provide also the theta parameter,\n";
    std::cerr << "  as well as min and max coord of the universe. example: \n";
    std::cerr << "  " << argv[0] << " [inputFilename] BarnesHut 0.5 -200 200\n";
    std::cerr << "NOTE! if using MPI solvers, run with: mpirun -n [procCount] ...\n";

    std::cerr << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printHelp(argv);
        return 1;
    }

    const std::string inputFilename = argv[1];
    const std::string solverType = argv[2];

    #if USE_OPENMP
        // Print OpenMP info
        #pragma omp parallel
        {
            #pragma omp master
            {
                std::cout << "Using OpenMP with " << omp_get_num_threads() << " threads" << std::endl;
            }
        }
    #endif


    // NOTE: Choose force.
    const forces::force<DIM, EmptyAttributes> &force = forces::gravity<DIM>();

    std::unique_ptr<NBodySolver<DIM, EmptyAttributes>> solver;
    
    if (solverType == "Serial") {
        solver = std::make_unique<SerialSolver<DIM, EmptyAttributes>>(force);
    } else if (solverType == "SerialReduced") {
        solver = std::make_unique<SerialReducedSolver<DIM, EmptyAttributes>>(force);
    } else if (solverType == "MPI") {
        solver = std::make_unique<MPISolver<DIM, EmptyAttributes>>(force);
    } else if (solverType == "MPIReduced") {
        solver = std::make_unique<MPIReducedSolver<DIM, EmptyAttributes>>(force);
    } else if (solverType == "BarnesHut") {
        if (argc < 6) {
            printHelp(argv);
            return 1;
        }
        const double bhTheta = std::stod(argv[3]);
        const double universeMinCoord = std::stod(argv[4]);
        const double universeMaxCoord = std::stod(argv[5]);
        const Box<DIM> universeBounds {Vec<DIM>(universeMinCoord), Vec<DIM>(universeMaxCoord)};
        solver = std::make_unique<BarnesHutSolver<DIM, EmptyAttributes>>(force, bhTheta, universeBounds);
    } else {
        std::cerr << "ERROR: UNKNOWN SOLVER TYPE " << solverType << std::endl;
        printHelp(argv);
        return 1;
    }

    // Init solver and accumulator.
    solver->init(argc, argv, inputFilename);
    AccelerationAccumulator<DIM, EmptyAttributes> &accumulator = solver->getAccumulator();

    // NOTE: Choose integrator.
    // integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
    integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

    // NOTE: Choose if you want output per timestep for visualization or not.
    // solver->runSimulation(integrator, solverType); // <--- use if no output steps (no vis!)
    solver->runSimulation(integrator, solverType, 1); // <--- use if output steps (then run vis tool)
    
}