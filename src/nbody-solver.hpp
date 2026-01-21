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

    void init(const std::string &filename) {
        readFromFile(filename);
    }

    void runSimulation(
        integrators::Integrator<DIM, Attributes> &integrator,
        bool printTimer=true
    ) {
        Timer t;
        t.start();

        for (int step = 0; step < steps; step++) {
            integrator.step(bodies, dt);
            oneStep();
        }

        t.end();
        if (printTimer)
            t.print();
    }

protected:
    int n; // body count
    int steps;
    double dt;
    Bodies<DIM, EmptyAttributes> bodies;
    std::unique_ptr<AccelerationAccumulator<DIM, Attributes>> accumulator;
    const forces::force<DIM, Attributes> &force;

    void oneStep() {} // Do something after integrator step if needed.

    void readFromFile(const std::string &filename) {
        utils::readFromFile<DIM>(filename, steps, dt, bodies);
        n = bodies.globalSize();
    }

};

template class NBodySolver<2, EmptyAttributes>;
template class NBodySolver<3, EmptyAttributes>;