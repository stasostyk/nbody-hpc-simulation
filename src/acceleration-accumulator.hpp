#pragma once

#include "body.hpp"

// Given a collection of bodies this class computes their accelerations.
template<int DIM>
class AccelerationAccumulator {
public:
    virtual void compute(bodies<DIM> &bodies) = 0;
    virtual const Vec<DIM>& accel(int bodyIndex) const = 0;
};
