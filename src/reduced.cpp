#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

constexpr int MAX_DIM = 3;
constexpr double G = 1.0;
constexpr double EPS = 1e-3;

struct Particle {
    double mass;
    double pos[MAX_DIM];
    double vel[MAX_DIM];
};

void compute_forces_reduced(const std::vector<Particle>& p,
                            std::vector<std::vector<double>>& forces,
                            int DIM)
{
    const int n = p.size();

    // azzera forze
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < DIM; ++d)
            forces[i][d] = 0.0;

    // solo coppie i < j
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {

            double dist2 = EPS * EPS;
            double dvec[MAX_DIM];

            for (int d = 0; d < DIM; ++d) {
                dvec[d] = p[i].pos[d] - p[j].pos[d];
                dist2 += dvec[d] * dvec[d];
            }

            double dist  = std::sqrt(dist2);
            double dist3 = dist2 * dist;

            double factor = G * p[i].mass * p[j].mass / dist3;

            // forza = factor * (ri - rj)
            for (int d = 0; d < DIM; ++d) {
                double f = factor * dvec[d];
                forces[i][d] -= f;
                forces[j][d] += f;  
            }
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


// EULER EXPLIXIT, update
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

int main()
{
    int DIM, n, n_steps;
    double dt;

    std::cin >> DIM >> n >> n_steps >> dt;

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

    auto tStart = std::chrono::steady_clock::now();


    // time loop
    for (int step = 0; step < n_steps; ++step) {

        // auto t0 = std::chrono::steady_clock::now();
        compute_forces_reduced(particles, forces, DIM);
        // auto t1 = std::chrono::steady_clock::now();

        // if (step == 0) {
        //     double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        //     std::cout << "Force compute time (ms): " << ms << "\n";
        // }

        // if (step % 200 == 0) {
        // std::cout << "Energy step " << step << ": " << total_energy(particles, DIM) << "\n";
        // }


        update_particles(particles, forces, DIM, dt);
    }

    auto tEnd = std::chrono::steady_clock::now();


    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart);

    std::cout << "Duration in ms: " << duration2.count() << std::endl;



    // // output finale
    // for (int i = 0; i < n; ++i) {
    //     std::cout << i << " " << particles[i].mass << " ";

    //     for (int d = 0; d < DIM; ++d)
    //         std::cout << particles[i].pos[d] << " ";

    //     for (int d = 0; d < DIM; ++d)
    //         std::cout << particles[i].vel[d] << " ";

    //     std::cout << "\n";
    // }

    return 0;
}
