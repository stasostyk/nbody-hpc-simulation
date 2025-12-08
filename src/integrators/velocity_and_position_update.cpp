#include "integrator.hpp"
#include "../forces_serial.hpp"

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
// Runge-Kutta 4 integrator (RK4) is much more accurate than Euler and Verlet
// it estimates the slope (derivative) four times per timestep
// and then takes a weighted average of those slopes
// this gives very good accuracy for short and medium time simulations
// but it's not symplectic, so energy can still slowly drift in long runs
void rk4(std::vector<Body> &bodies, double dt){
    size_t n = bodies.size();
    // k?_s are the slopes for position (ds/dt = v)
    // k?_v are the slopes for velocity (dv/dt = a)
    std::vector<vect_t> k1_s(n), k2_s(n), k3_s(n), k4_s(n);
    std::vector<vect_t> k1_v(n), k2_v(n), k3_v(n), k4_v(n);
    std::vector<Body> temp = bodies; // Temporary copy for intermediate steps
    
    // Step 1: Compute k1
    compute_acceleration(temp);
    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            k1_s[i][d] = temp[i].v[d];
            k1_v[i][d] = temp[i].a[d];
        }
    }

    // Step 2: Compute k2
    // we move the system forward by dt/2 using k1 slopes
    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            temp[i].s[d] = bodies[i].s[d] + (0.5 * dt) * k1_s[i][d];
            temp[i].v[d] = bodies[i].v[d] + (0.5 * dt) * k1_v[i][d];
        }
    }

    compute_acceleration(temp);  // a at s^(2)

    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            k2_s[i][d] = temp[i].v[d];
            k2_v[i][d] = temp[i].a[d];
        }
    }

    // Step 3: compute k3, again at dt/2 but now using k2 slopes
    // this improves the midpoint estimate
    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            temp[i].s[d] = bodies[i].s[d] + (0.5 * dt) * k2_s[i][d];
            temp[i].v[d] = bodies[i].v[d] + (0.5 * dt) * k2_v[i][d];
        }
    }

    compute_acceleration(temp);

    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            k3_s[i][d] = temp[i].v[d];
            k3_v[i][d] = temp[i].a[d];
        }
    }
    // Step 4: compute k4 using a full timestep with k3
    // this gives the slope at the end of the timestep
    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            temp[i].s[d] = bodies[i].s[d] + dt * k3_s[i][d];
            temp[i].v[d] = bodies[i].v[d] + dt * k3_v[i][d];
        }
    }

    compute_acceleration(temp);

    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            k4_s[i][d] = temp[i].v[d];
            k4_v[i][d] = temp[i].a[d];
        }
    }
    
    // Final and actual update
    // we combine all four slopes using RK4 weights:
    // k1 + 2*k2 + 2*k3 + k4, divided by 6
    for (size_t i = 0; i < n; ++i) {
        for (int d = 0; d < DIM; ++d) {
            bodies[i].v[d] += (dt/6.0) * (k1_v[i][d] + 2.0 * k2_v[i][d] + 2.0 * k3_v[i][d] + k4_v[i][d]);
            bodies[i].s[d] += (dt/6.0) * (k1_s[i][d] + 2.0 * k2_s[i][d] + 2.0 * k3_s[i][d] + k4_s[i][d]);
        }
    }
}
