#include "BHTree.hpp"
#include "../utils.hpp"
#include "../forces/gravity.hpp"
#include "../timer.hpp"
#include "bhAccumulator.hpp"
#include "../integrators/euler.hpp"
#include "../integrators/integrator.hpp"
#include "../integrators/sympletic.hpp"
#include "../integrators/verlet.hpp"
#include "../integrators/rk4.hpp"
#include "../acceleration-accumulator.hpp"

#if USE_OPENMP
    #include <omp.h>
#endif

const int DIM = 2;

int main(int argc, char* argv[]) {

    // utils::compareOutputs<DIM>("serial.test1.out", "bhSerial.test1.lastStep.out");
    // return 0;

    double theta = 0.5; // theta parameter for BH

    if (argc > 1) {
        theta = std::stod(argv[1]);
    }

    Box<DIM> universeBounds(Vec<DIM>(-200.), Vec<DIM>(200.));

    const forces::force<DIM, EmptyAttributes> &force = forces::gravity<DIM>();

    int n;
    int steps;
    double dt;
    bodies<DIM, EmptyAttributes> bodies;

    utils::generateRandomToFile<DIM>("test1.in.out", 1000, 1000, 0.01, 42);
    utils::readFromFile("test1.in.out", steps, dt, bodies);
    n = bodies.globalSize();

    BarnesHutAccumulator<DIM, EmptyAttributes> accumulator(n, theta, universeBounds, force);

    integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);


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

    Timer t;
    t.start();

    for (int step = 0; step < steps; step++) {
        integrator.step(bodies, dt);
    }

    t.end();
    t.print();

    // save for testing
    utils::saveToFile<DIM>("bhSerial.test1.lastStep.out", steps, dt, bodies, false);
}