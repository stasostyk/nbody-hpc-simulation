#include <memory>

#include "nbody-solver.hpp"
#include "serialSolvers/serialSolver.hpp"
#include "serialSolvers/serialReducedSolver.hpp"
#include "mpiSolvers/mpiSolver.hpp"
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

#define DIM 2

void printHelp(char **argv) {
    std::cerr << "Usage: " << argv[0]
                << " [inputFilename] [solver]\n";
    std::cerr << "  solver: Serial SerialReduced MPI\n";
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


    // NOTE: Choose force.
    const forces::force<DIM, EmptyAttributes> &force = forces::gravity<DIM>();

    std::unique_ptr<NBodySolver<DIM, EmptyAttributes>> solver;
    
    if (solverType == "Serial") {
        solver = std::make_unique<SerialSolver<DIM, EmptyAttributes>>(force);
    } else if (solverType == "SerialReduced") {
        solver = std::make_unique<SerialReducedSolver<DIM, EmptyAttributes>>(force);
    } else if (solverType == "MPI") {
        solver = std::make_unique<MPISolver<DIM, EmptyAttributes>>(force);
    } else {
        std::cerr << "ERROR: UNKNOWN SOLVER TYPE " << solverType << std::endl;
        printHelp(argv);
        return 1;
    }

    // SerialSolver<DIM, EmptyAttributes> solver(force);
    
    // Init solver and accumulator.
    solver->init(argc, argv, inputFilename);
    AccelerationAccumulator<DIM, EmptyAttributes> &accumulator = solver->getAccumulator();

    // NOTE: Choose integrator.
    // TODO choose integrator via command line param
    integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

    solver->runSimulation(integrator, solverType);
}