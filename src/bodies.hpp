#pragma once
#include <array>
#include "vec.hpp"

template<size_t dim>
struct Body
{
    unsigned int bodyId; // to distinguish bodies, needed in BHTree code

    double mass;
    Vec<dim> position;
    Vec<dim> velocity;
    Vec<dim> acceleration;
    Vec<dim> force;

    void printPosition() const {
        std::cout << "(";
        for (size_t dimId = 0; dimId < dim-1; dimId++)
            std::cout << position[dimId] << ", ";
        std::cout << position[dim-1] << ")";
    }
};
