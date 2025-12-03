#include "integrator.hpp"
#include "forces_serial.hpp"
#include <vector>
#include "bodies.hpp"
// basic euler integrator it's the simplest but it can drift from accuracy easily
//and overshhot or undershhot aproximations
//it works by updating the position using OLD velocity (wgy it's not accurate)
//then updating velocity using acceleration
void euler(std::vector<Body> &bodies, double dt){
    compute_acceleration(bodies);
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


// sympletic euler integrator it's more accurate than basic euler
//diffrience is that it updates velocity first then position using NEW velocity
//that's why it's more accurate
void sympletic_euler(std::vector<Body> &bodies, double dt){
    compute_acceleration(bodies);
    size_t n = bodies.size();
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            bodies[i].v[d] += bodies[i].a[d] * dt;
            bodies[i].s[d] += bodies[i].v[d] * dt;
        }
    }
}

