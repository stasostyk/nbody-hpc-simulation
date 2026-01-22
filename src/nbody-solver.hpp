#pragma once

#include <memory>

#include "forces/func.hpp"
#include "integrators/integrator.hpp"
#include "timer.hpp"
#include "utils.hpp"
#include "acceleration-accumulator.hpp"

template <int DIM, typename Attributes>
class NBodySolver {
public:
    NBodySolver(
        const forces::force<DIM, Attributes> &_force
    ) : force(_force) {}

    virtual AccelerationAccumulator<DIM, Attributes> &getAccumulator() = 0;

    virtual void init(int argc, char **argv, const std::string &filename) {
        readFromFile(filename);
    }

    virtual void finalize() {} // in MPI Solver case, they need to finalize

    virtual void runSimulation(
        integrators::Integrator<DIM, Attributes> &integrator,
        const std::string &outputFilePref,
        int outputFrequency = -1
    ) {
        Timer t;
        t.start();

        for (int step = 0; step < steps; step++) {
            integrator.step(bodies, dt);
            oneStep();

            if (outputFrequency != -1 && (step % outputFrequency) == 0) {
                std::string outputFinalFilename = 
                    outputFilePref + ".step-" + std::to_string(step / outputFrequency) + ".out";
                saveOutput(outputFinalFilename);
            }
        }

        t.end();
        t.print();

        std::string outputFinalFilename = outputFilePref + ".finalTimestep.out";
        saveOutput(outputFinalFilename);

        finalize();
    }

protected:
    int n; // body count
    int steps;
    double dt;
    Bodies<DIM, EmptyAttributes> bodies;
    std::unique_ptr<AccelerationAccumulator<DIM, Attributes>> accumulator;
    const forces::force<DIM, Attributes> &force;

    virtual void oneStep() {} // Do something after integrator step if needed.

    void readFromFile(const std::string &filename) {
        utils::readFromFile<DIM>(filename, steps, dt, bodies);
        n = bodies.globalSize();
    }

    virtual void saveOutput(const std::string &filename) {
        utils::saveToFile<DIM>(filename, steps, dt, bodies, false);
    }
    
};
