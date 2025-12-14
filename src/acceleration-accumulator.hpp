#pragma once

#include "body.hpp"

template<int DIM>
class AccelerationAccumulator {
public:
    virtual void compute(const std::vector<body<DIM>> &bodies) const = 0;
    virtual const Vec<DIM>& accel(int particleIndex) const = 0;
};
