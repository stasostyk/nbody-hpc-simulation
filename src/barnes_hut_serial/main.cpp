#include "BHTree.hpp"
#include "../utils.hpp"
#include "../forces/gravity.hpp"
#include "../timer.hpp"

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

    forces::gravity<DIM> gravity{};
    const forces::force<DIM> &force = gravity;

    int n;
    int steps;
    double dt;

    std::vector<double> masses;
    std::vector<Vec<DIM>> positions;
    std::vector<Vec<DIM>> velocities;
    std::vector<Vec<DIM>> forces;

    utils::generateRandomToFile<DIM>("test1.in.out", 1000, 1000, 0.01, 42);
    utils::readFromFile<DIM>("test1.in.out", n, steps, dt, masses, positions, velocities);

    std::vector<Body<DIM>> bodies(n);
    for (int i = 0; i < n; i++) {
        bodies[i].bodyId = i;
        bodies[i].acceleration = 0.;
        bodies[i].force = 0.;
        bodies[i].mass = masses[i];
        bodies[i].position = positions[i];
        bodies[i].velocity = velocities[i];
    }

    forces.resize(n);

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
        BHTree<DIM> bhTree(universeBounds, bodies);

        // Compute all forces using BH Tree.
        #if USE_OPENMP
            #pragma omp parallel for schedule(dynamic, 16)
        #endif
        for (int i = 0; i < n; i++) {
            forces[i] = bhTree.calculateForce(bodies[i], theta, force);
        }
 
        // apply forces, i.e. calc new positions and velocities
        #if USE_OPENMP
            #pragma omp parallel for schedule(static)
        #endif
        for (int i = 0; i < n; i++) {
            positions[i] += dt * velocities[i];
            velocities[i] += dt / masses[i] * forces[i];
        }
    
        // // save for testing
        // utils::saveToFile("test1." + std::to_string(step) + ".out", n, steps, dt,
        //     masses, positions, velocities, false);
        // if (step % 100 == 0) std::cout << " step " << step << " finished" << std::endl;
    }

    t.end();
    t.print();

    // save for testing
    utils::saveToFile("bhSerial.test1.lastStep.out", n, steps, dt,
        masses, positions, velocities, false);


    // // SOME OLD TESTS BELOW
    // std::vector<Body<2>> bodies(3);
    // bodies[0].position = Vec<2>({0.0, 0.0});
    // bodies[1].position = Vec<2>({0.5, 1.1});
    // bodies[2].position = Vec<2>({0.5, 1.4});

    // bodies[0].mass = 1.;
    // bodies[1].mass = 2.;
    // bodies[2].mass = 3.;

    // const Box<2> universeBounds(Vec<2>({0.0, 0.0}), Vec<2>({2.0, 2.0}));
    // BHTree<2> bhTree(universeBounds, bodies);

    // bhTree.printInfo();

    // Box<2> tmpBox (Vec<2>({0., 0.}), Vec<2>({1., 1.}));
    // std::cout << "is inside? " << tmpBox.isPointInside(Vec<2>(bodies[0].position)) << std::endl;
}