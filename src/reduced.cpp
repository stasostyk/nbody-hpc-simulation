#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>  

constexpr int DIM = 2;
 
#include "forces/gravity.hpp"

struct Particle {
    double mass;
    Vec<DIM> pos;
    Vec<DIM> vel;
};

void compute_forces_reduced(const std::vector<Particle>& p,
                            std::vector<Vec<DIM>>& forces, const forces::force<DIM>& force)
{
    int n = p.size();

    // azzera forze
    for (int q = 0; q < n; ++q) {
        forces[q][0] = 0.0;
        forces[q][1] = 0.0;
    }

    for (int q = 0; q < n; ++q) {
        for (int k = q + 1; k < n; ++k) {
            auto f = force(p[q].pos, p[q].mass, p[k].pos, p[k].mass);
            forces[q] += f;
            forces[k] -= f;
        }
    }
}


/* SEMBRA CHE QUESTA ABBIA QUALCHE SEGNO SBAGLIATO

void compute_forces_reduced(
    const std::vector<Particle>& p,
    std::vector<vect_t>& forces)
{
    int n = p.size();

    // reset
    for (int q = 0; q < n; ++q) {
        forces[q][0] = 0.0;
        forces[q][1] = 0.0;
    }

    for (int q = 0; q < n; ++q) {
        for (int k = q + 1; k < n; ++k) {

            double dx = p[q].pos[0] - p[k].pos[0];
            double dy = p[q].pos[1] - p[k].pos[1];

            double eps = 1e-9; // softening
            double dist2 = dx*dx + dy*dy + eps*eps;
            double dist  = std::sqrt(dist2);
            double dist3 = dist2 * dist;


            double factor = - G * p[q].mass * p[k].mass / dist3; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            double fx = factor * dx;
            double fy = factor * dy;

            forces[q][0] += fx;
            forces[q][1] += fy;

            forces[k][0] -= fx;
            forces[k][1] -= fy;
        }
    }
}*/

/* IN TEORIA QUESTO EVITA I NaN

void compute_forces_reduced(const std::vector<Particle>& p,
                            std::vector<vect_t>& forces)
{
    int n = p.size();

    for (int q = 0; q < n; ++q) {
        forces[q][0] = 0.0;
        forces[q][1] = 0.0;
    }

    const double eps = 1e-9;

    for (int q = 0; q < n; ++q) {
        for (int k = q + 1; k < n; ++k) {

            double dx = p[q].pos[0] - p[k].pos[0];
            double dy = p[q].pos[1] - p[k].pos[1];

            double dist2 = dx*dx + dy*dy + eps*eps;
            double dist  = std::sqrt(dist2);
            double dist3 = dist2 * dist;

            double factor = G * p[q].mass * p[k].mass / dist3;

            double fx = factor * dx;
            double fy = factor * dy;

            forces[q][0] += fx;
            forces[q][1] += fy;

            forces[k][0] -= fx;
            forces[k][1] -= fy;
        }
    }
}*/


void update_particles(std::vector<Particle>& p,
                      const std::vector<Vec<DIM>>& forces,
                      double dt)
{
    int n = p.size();
    int step = 0;
    for (int q = 0; q < n; ++q) {

        p[q].pos[0] += dt * p[q].vel[0];
        p[q].pos[1] += dt * p[q].vel[1];

        p[q].vel[0] += dt * forces[q][0] / p[q].mass;
        p[q].vel[1] += dt * forces[q][1] / p[q].mass;

        if (std::isnan(p[q].pos[0]) || std::isnan(p[q].pos[1]) ||
            std::isnan(p[q].vel[0]) || std::isnan(p[q].vel[1])) {

            std::cerr << "NaN at particle " << q << " step: " << step << "\n";
            exit(1);
        }
        step++;

    }
}

int main(int argc, char* argv[]) {

    //      STRESS TEST 4 (many particles)

    /*if (argc >= 2 && strcmp(argv[1], "--test4") == 0) {

        // default N = 200
        int N = 200;
        if (argc >= 3) {
            N = std::atoi(argv[2]);
            if (N <= 0) N = 200;
        }

        int n = N;
        int n_steps = 1000;
        double dt = 0.001;

        std::cout << n << " " << n_steps << " " << dt << "\n";

        std::mt19937 rng(123);
        std::uniform_real_distribution<double> Upos(-1.0, 1.0);
        std::uniform_real_distribution<double> Uvel(-0.5, 0.5);

        for (int i = 0; i < n; ++i) {
            double mass = 1.0;
            double x = Upos(rng);
            double y = Upos(rng);
            double vx = Uvel(rng);
            double vy = Uvel(rng);

            std::cout << mass << " "
                      << x << " " << y << " "
                      << vx << " " << vy << "\n";
        }

        return 0;  
    }*/

    int dim, n, n_steps;
    double dt;
    std::cin >> dim >> n >> n_steps >> dt;

    std::vector<Particle> particles(n);
    std::vector<Vec<DIM>> forces(n);

    for (int q = 0; q < n; ++q) {
        std::cin >> particles[q].mass
                 >> particles[q].pos[0] >> particles[q].pos[1]
                 >> particles[q].vel[0] >> particles[q].vel[1];
    }

    for (int step = 0; step < n_steps; ++step) {
        compute_forces_reduced(particles, forces, forces::gravity<DIM>{});
        update_particles(particles, forces, dt);
    }

    for (int q = 0; q < n; ++q) {
        std::cout << q << " "
                  << particles[q].mass << " "
                  << particles[q].pos[0] << " " << particles[q].pos[1] << " "
                  << particles[q].vel[0] << " " << particles[q].vel[1] << "\n";
    }

    return 0;
}
