#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>


constexpr int DIM = 2;
constexpr double G = 1.0;
constexpr double EPS = 1e-3; //1e-6

using vect_t = double[DIM];

struct Particle {
    double mass;
    vect_t pos;
    vect_t vel;
};

inline void compute_forces_basic(const std::vector<Particle>& p,
                                 std::vector<vect_t>& forces)
{
    const int n = p.size();

    // Azzera forze
    for (int q = 0; q < n; ++q) {
        forces[q][0] = 0.0;
        forces[q][1] = 0.0;
    }

    // Tutte le coppie q != k
    for (int q = 0; q < n; ++q) {
        for (int k = 0; k < n; ++k) {
            if (q == k) continue;

            const double dx = p[q].pos[0] - p[k].pos[0];
            const double dy = p[q].pos[1] - p[k].pos[1];

            const double dist2 = dx*dx + dy*dy + EPS*EPS;
            const double dist  = std::sqrt(dist2);
            const double dist3 = dist2 * dist;

            const double factor = G * p[q].mass * p[k].mass / dist3;

            if (std::isnan(dist2) || std::isnan(dist3)) {
            std::cout << "NaN at q=" << q << ", k=" << k 
              << " dx=" << dx << " dy=" << dy << "\n";
            exit(1);
            }

            forces[q][0] -= factor * dx;
            forces[q][1] -= factor * dy;
        }
    }
}

inline void update_particles(std::vector<Particle>& p,
                             const std::vector<vect_t>& forces,
                             double dt)
{
    const int n = p.size();

    for (int q = 0; q < n; ++q) {

        // Update posizione (x(t+dt) = x(t) + dt*v(t))
        p[q].pos[0] += dt * p[q].vel[0];
        p[q].pos[1] += dt * p[q].vel[1];

        // Update velocità (v(t+dt) = v(t) + dt*a)
        p[q].vel[0] += dt * forces[q][0] / p[q].mass;
        p[q].vel[1] += dt * forces[q][1] / p[q].mass;
    }
}

//=========================== VERIFICA CONSERVAZIONE ENERGIA ==============================
double total_energy(const std::vector<Particle>& p) {
    double E = 0.0;
    int n = p.size();

    // Energia cinetica
    for (int i = 0; i < n; i++) {
        double vx = p[i].vel[0];
        double vy = p[i].vel[1];
        E += 0.5 * p[i].mass * (vx*vx + vy*vy);
    }

    // Energia potenziale gravitazionale
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            double dx = p[i].pos[0] - p[j].pos[0];
            double dy = p[i].pos[1] - p[j].pos[1];
            double dist = std::sqrt(dx*dx + dy*dy + 1e-6*1e-6);

            E -= G * p[i].mass * p[j].mass / dist;
        }
    }

    return E;
}


int main()
{
    int n, n_steps;
    double dt;

    std::cin >> n >> n_steps >> dt;

    std::vector<Particle> particles(n);
    std::vector<vect_t> forces(n);

    int k=1;

    // Input particelle
    for (int q = 0; q < n; ++q) {
        std::cin >> particles[q].mass
                 >> particles[q].pos[0] >> particles[q].pos[1]
                 >> particles[q].vel[0] >> particles[q].vel[1];
        k++;
    }

    if(k!=n){
        std::cout << "MAKE SURE THAT N=NUMEBER_OF_PARTICLES\n\n";
    }


    for (int q = 0; q < n; ++q) {
    if (std::isnan(particles[q].pos[0]) ||
        std::isnan(particles[q].pos[1]) ||
        std::isnan(particles[q].vel[0]) ||
        std::isnan(particles[q].vel[1])) 
    {
        std::cout << "NaN in input at particle " << q << "\n";
        exit(1);
    }
    }


    // Time loop
    for (int step = 0; step < n_steps; ++step) {
        

            // --- TIMER START ---
        auto t0 = std::chrono::steady_clock::now();
        compute_forces_basic(particles, forces);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        // --- TIMER END ---



        if (step == 0) std::cout << "Force compute time (ms): " << ms << "\n";

        update_particles(particles, forces, dt);
        /*if (step % 1000 == 0) {
            std::cout << "Energy step " << step << ": " << total_energy(particles) << "\n";
        }*/

    }

    // Output finale
    for (int q = 0; q < n; ++q) {
        std::cout << q << " "
                  << particles[q].mass << " "
                  << particles[q].pos[0] << " " << particles[q].pos[1] << " "
                  << particles[q].vel[0] << " " << particles[q].vel[1] << "\n";
    }

    return 0;
}
