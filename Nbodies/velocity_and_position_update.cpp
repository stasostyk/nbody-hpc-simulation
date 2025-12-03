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

// velocity verlet integrator is even more accurate
// it updates velocity in half steps and then use that half step v to compute position update
//recomputes acceleration using updated positions
// this way it takes into account the change in acceleration due to position change
void velocity_verlet(std::vector<Body> &bodies, double dt){
    compute_acceleration(bodies);
    size_t n = bodies.size();
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            bodies[i].v[d] += bodies[i].a[d] * 0.5 * dt; //we update velocity only by half step
            bodies[i].s[d] += bodies[i].v[d] * dt; //then we update position using this half step velocity
        }
    }
    compute_acceleration(bodies);       //recompute acceleration with new positions
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            bodies[i].v[d] += bodies[i].a[d] * 0.5 * dt;    //complete velocity update with the new acceleration
        }
    }
}
