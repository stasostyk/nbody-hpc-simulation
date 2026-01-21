#pragma once

#include "../nbody-solver.hpp"
#include "serial-accumulator-reduced.hpp"

template <int DIM, typename Attributes>
class SerialReducedSolver : public NBodySolver<DIM, Attributes>
{
public:
    using NBodySolver<DIM, Attributes>::NBodySolver;  // Inherit constructors

    AccelerationAccumulator<DIM, Attributes> &getAccumulator() override
    {
        this->accumulator = std::make_unique<SerialAccumulatorReduced<DIM, Attributes>>(
            this->bodies.localSize(), this->force
        );
        return *(this->accumulator);
    }
};
