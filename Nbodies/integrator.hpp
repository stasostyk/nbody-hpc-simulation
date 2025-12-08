#pragma once
#include <vector>
#include "bodies.hpp"
void euler(std::vector<Body> &bodies, double dt);
void sympletic_euler(std::vector<Body> &bodies, double dt);
void velocity_verlet(std::vector<Body> &bodies, double dt);
void rk4(std::vector<Body> &bodies, double dt);