#include <memory>

#include "nbody-solver.hpp"
#include "serialSolvers/serialSolver.hpp"
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

int main() {
    const std::string inputFilename = "test1.in.out";
    utils::generateRandomToFile<DIM>(inputFilename, 1000, 1000, 0.01, 42);

    // NOTE: Choose force.
    const forces::force<DIM, EmptyAttributes> &force = forces::gravity<DIM>();

    // NOTE: Choose solver.
    SerialSolver<DIM, EmptyAttributes> solver(force);
    
    // Init solver and accumulator.
    solver.init(inputFilename);
    AccelerationAccumulator<DIM, EmptyAttributes> &accumulator = solver.getAccumulator();

    // NOTE: Choose integrator.
    integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);

    solver.runSimulation(integrator, true);
}