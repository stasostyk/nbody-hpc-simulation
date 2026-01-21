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

    void combineWith(const Body<dim> &other) {
        if (other.mass < 1e-10) return;

        double newMass = mass + other.mass;
        Vec<dim> newPosition = (position * mass + other.position * other.mass) / newMass;
        Vec<dim> newVelocity = (velocity * mass + other.velocity * other.mass) / newMass;

        mass = newMass;
        position = newPosition;
        velocity = newVelocity;
    }
};
