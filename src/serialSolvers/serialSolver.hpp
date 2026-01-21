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
        // return SerialAccumulator<DIM, Attributes>(this->bodies.localSize(), this->force);
    }
};

template class SerialSolver<2, EmptyAttributes>;
template class SerialSolver<3, EmptyAttributes>;
