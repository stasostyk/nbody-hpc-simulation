#pragma once

#include "../nbody-solver.hpp"
#include "serial-accumulator.hpp"

template <int DIM, typename Attributes>
class SerialSolver : public NBodySolver<DIM, Attributes>
{
public:
    using NBodySolver<DIM, Attributes>::NBodySolver;  // Inherit constructors

    AccelerationAccumulator<DIM, Attributes> &getAccumulator() override
    {
        this->accumulator = std::make_unique<SerialAccumulator<DIM, Attributes>>(
            this->bodies.localSize(), this->force
        );
        return *(this->accumulator);
    }
};
