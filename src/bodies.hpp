#pragma once

#include <array>

#define DIM 2
using vect_t = std::array<double, DIM>;
struct Body
{
    double mass;
    vect_t s; // position
    vect_t v; // velocity
    vect_t a; // acceleration
    vect_t f; // force
};
