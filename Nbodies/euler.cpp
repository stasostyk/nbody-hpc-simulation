#include "integrator.hpp"
#include "forces.hpp"
#include <vector>
#include "bodies.hpp"
void euler(std::vector<Body> &bodies, double dt){
    compute_acceleration(bodies); // compute acceleration at time t
    size_t n = bodies.size();
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            bodies[i].s[d] += bodies[i].v[d] * dt;
        }
    }
        for (int i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            bodies[i].v[d] += bodies[i].a[d] * dt;
        }
    }
}