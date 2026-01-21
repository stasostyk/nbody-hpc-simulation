#include "../nbody-solver.hpp"
#include "bhAccumulator.hpp"

template <int DIM, typename Attributes>
class BarnesHutSolver : public NBodySolver<DIM, Attributes>
{
public:
    using NBodySolver<DIM, Attributes>::NBodySolver;  // Inherit constructors
    
    BarnesHutSolver(
        const forces::force<DIM, Attributes> &_force,
        const double _theta,
        const Box<DIM> &_universeBounds
    ) : NBodySolver<DIM, Attributes>(_force),
        theta(_theta), universeBounds(_universeBounds) {}

    AccelerationAccumulator<DIM, Attributes> &getAccumulator() override
    {
        this->accumulator = std::make_unique<BarnesHutAccumulator<DIM, Attributes>>(
            this->n, theta, universeBounds, this->force
        );
        return *(this->accumulator);
    }
private:
    const double theta;
    const Box<DIM> universeBounds;
};
