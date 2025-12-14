#pragma once

#include "vec.hpp"

template<int DIM>
struct body
{
    Vec<DIM> position;
    Vec<DIM> velocity;
    double mass;
};
