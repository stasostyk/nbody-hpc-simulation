#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>

constexpr int MAX_DIM = 3;
constexpr double G = 1.0;
constexpr double EPS = 1e-3;

struct Particle {
    double mass;
    double pos[MAX_DIM];
    double vel[MAX_DIM];
};

void compute_forces_basic(const std::vector<Particle>& p,
                          std::vector<std::vector<double>>& forces,
                          int DIM)
{
    const int n = p.size();

    // azzera forze
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < DIM; ++d)
            forces[i][d] = 0.0;

    // tutte le coppie q != k
    for (int q = 0; q < n; ++q) {
        for (int k = 0; k < n; ++k) {
            if (q == k) continue;

            double dist2 = EPS * EPS;
            double dvec[MAX_DIM];

            for (int d = 0; d < DIM; ++d) {
                dvec[d] = p[q].pos[d] - p[k].pos[d];
                dist2 += dvec[d] * dvec[d];
            }

            double dist  = std::sqrt(dist2);
            double dist3 = dist2 * dist;

            double factor = G * p[q].mass * p[k].mass / dist3;

            for (int d = 0; d < DIM; ++d)
                forces[q][d] -= factor * dvec[d];
        }
    }
}

double total_energy(const std::vector<Particle>& p, int DIM) {
    const int n = (int)p.size();
    double E = 0.0;

    // Kinetic
    for (int i = 0; i < n; ++i) {
        double v2 = 0.0;
        for (int d = 0; d < DIM; ++d) v2 += p[i].vel[d] * p[i].vel[d];
        E += 0.5 * p[i].mass * v2;
    }

    // Potential (i<j)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double r2 = 0.0;
            for (int d = 0; d < DIM; ++d) {
                double dx = p[i].pos[d] - p[j].pos[d];
                r2 += dx * dx;
            }
            double r = std::sqrt(r2 + EPS*EPS);
            E -= G * p[i].mass * p[j].mass / r;
        }
    }
    return E;
}

// EULER EXPLICIT
void update_particles(std::vector<Particle>& p,
                      const std::vector<std::vector<double>>& forces,
                      int DIM, double dt)
{
    const int n = p.size();

    for (int i = 0; i < n; ++i) {

        // posizione
        for (int d = 0; d < DIM; ++d)
            p[i].pos[d] += dt * p[i].vel[d];

        // velocità
        for (int d = 0; d < DIM; ++d)
            p[i].vel[d] += dt * forces[i][d] / p[i].mass;
    }
}

double critical_timestep(const Particle& p,
                         const std::vector<double>& force,
                         int DIM,
                         double dt_max) {
    double acc_norm_sq = 0.0;
    double vel_norm_sq = 0.0;
    for (int d = 0; d < DIM; ++d) {
        double acc = force[d] / p.mass;
        acc_norm_sq += acc * acc;
        vel_norm_sq += p.vel[d] * p.vel[d];
    }

    double acc_norm = std::sqrt(acc_norm_sq);
    double vel_norm = std::sqrt(vel_norm_sq);

    // Use the softening length as characteristic length scale for stability.
    double acc_limit = acc_norm > 0.0 ? std::sqrt(EPS / acc_norm) : dt_max;
    double vel_limit = vel_norm > 0.0 ? EPS / vel_norm : dt_max;

    double dt_crit = std::min({dt_max, acc_limit, vel_limit});
    const double min_dt = dt_max * 1e-6; // prevent vanishing timesteps
    return std::max(dt_crit, min_dt);
}

double update_timesteps(const std::vector<Particle>& particles,
                        const std::vector<std::vector<double>>& forces,
                        int DIM,
                        double dt_max,
                        std::vector<double>& particle_dt) {
    const double min_dt = dt_max * 1e-6;
    double next_dt = dt_max;

    for (size_t i = 0; i < particles.size(); ++i) {
        double dt_crit = critical_timestep(particles[i], forces[i], DIM, dt_max);
        double candidate = particle_dt[i];

        if (candidate > dt_crit) {
            while (candidate > dt_crit) {
                candidate *= 0.5; // refine until stable
            }
        } else if (candidate < 0.5 * dt_crit) {
            candidate = std::min(candidate * 2.0, dt_max); // coarsen
        }

        candidate = std::clamp(candidate, min_dt, dt_max);
        particle_dt[i] = candidate;
        next_dt = std::min(next_dt, candidate);
    }

    return next_dt;
}


int main()
{
    int DIM, n;
    double total_time, dt_max;

    std::cin >> DIM >> n >> total_time >> dt_max;

    if (DIM < 2 || DIM > 3) {
        std::cerr << "DIM must be 2 or 3\n";
        return 1;
    }

    std::vector<Particle> particles(n);
    std::vector<std::vector<double>> forces(n, std::vector<double>(DIM));

    // input particelle
    for (int i = 0; i < n; ++i) {
        std::cin >> particles[i].mass;

        for (int d = 0; d < DIM; ++d)
            std::cin >> particles[i].pos[d];

        for (int d = 0; d < DIM; ++d)
            std::cin >> particles[i].vel[d];
    }

    double time = 0.0;
    int step = 0;
    std::vector<double> particle_dt(n, dt_max);

    // time loop
    while (time < total_time) {

        auto t0 = std::chrono::steady_clock::now();
        compute_forces_basic(particles, forces, DIM);
        auto t1 = std::chrono::steady_clock::now();

        if (step == 0) {
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::cout << "Force compute time (ms): " << ms << "\n";
        }

        if (step % 200 == 0) {
            std::cout << "Energy t=" << time << " (step " << step << "): " << total_energy(particles, DIM) << "\n";
        }

        double dt = update_timesteps(particles, forces, DIM, dt_max, particle_dt);
        dt = std::min(dt, total_time - time);
        if (dt <= 0.0) break;

        update_particles(particles, forces, DIM, dt);

        time += dt;
        ++step;
    }

    // output finale
    for (int i = 0; i < n; ++i) {
        std::cout << i << " " << particles[i].mass << " ";

        for (int d = 0; d < DIM; ++d)
            std::cout << particles[i].pos[d] << " ";

        for (int d = 0; d < DIM; ++d)
            std::cout << particles[i].vel[d] << " ";

        std::cout << "\n";
    }

    return 0;
}
