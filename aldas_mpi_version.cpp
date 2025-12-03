/*

compilation:
mpic++ -O3 -std=c++11 -o nbody_mpi_reduced_packed nbody_mpi_reduced_packed.cpp

running:
mpirun -n 4 ./nbody_mpi_reduced_packed 400 100 0.01 0

*/

#include <assert.h>
#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>

const static int DIM = 2; // dimensions of the problem (2D or 3D)
const double G = 6.673e-11; // gravitational constant

using Vec = std::array<double, DIM>; // vector in our DIM-dimensional space

// struct Particle {
//     double mass;
//     double velocity;
//     std::array<double, DIM> position;
//     std::array<double, DIM> force;

//     Particle() : mass(0.), velocity(0.), position{}, force{} {}
// };

// // gets force acting on particle p1
// inline std::array<double, DIM> force_on_p1(const Particle &p1, const Particle &p2)
// {
//     std::array<double, DIM> dist;
//     double distSquared = 1e-12; // softening to avoid division by zero
//     for (int i = 0; i < DIM; i++) {
//         dist[i] = p1.position[i] - p2.position[i];
//         distSquared += dist[i] * dist[i];
//     }

//     double dist3 = distSquared * std::sqrt(distSquared);
//     double coeff = -G * p1.mass * p2.mass / dist3;

//     std::array<double, DIM> force;
//     for (int i = 0; i < DIM; i++) 
//         force[i] = dist[i] * coeff;
    
//     return force;
// }


// gets force acting on particle p1
inline Vec force_on_p1(
    const Vec &pos1, 
    const Vec &pos2,
    double m1, double m2
) {
    Vec dist;
    double distSquared = 1e-12; // softening to avoid division by zero
    for (int i = 0; i < DIM; i++) {
        dist[i] = pos1[i] - pos2[i];
        distSquared += dist[i] * dist[i];
    }

    double dist3 = distSquared * std::sqrt(distSquared);
    double coeff = -G * m1 * m2 / dist3;

    Vec force;
    for (int i = 0; i < DIM; i++) 
        force[i] = dist[i] * coeff;
    
    return force;
}

inline void saveToFile(
    const std::string &fileName,
    int n, int steps, double dt,
    const std::vector<double> &m, 
    const std::vector<Vec> &p, 
    const std::vector<Vec> &v
) {
    std::ofstream fout(fileName);

    fout << DIM << "\n";
    fout << n << " " << steps << " " << dt << "\n";
    for (int i = 0; i < n; i++) {
        fout << m[i] << " ";
        for (int j = 0; j < DIM; j++) fout << p[i][j] << " ";
        for (int j = 0; j < DIM; j++) fout << v[i][j] << " ";
        fout << "\n";
    }
    fout.flush();
    fout.close();

    std::cout << "Saved to file " << fileName << std::endl;
}

inline void readFromFile(
    const std::string &fileName,
    int &n, int &steps, double &dt,
    std::vector<double> &m, 
    std::vector<Vec> &p, 
    std::vector<Vec> &v
) {
    std::ifstream fin(fileName);

    int dimInFile;
    fin >> dimInFile;
    assert(dimInFile == DIM);

    fin >> n >> steps >> dt;

    m.resize(n);
    p.resize(n);
    v.resize(n);

    for (int i = 0; i < n; i++) {
        fin >> m[i];
        for (int j = 0; j < DIM; j++) fin >> p[i][j];
        for (int j = 0; j < DIM; j++) fin >> v[i][j];
    }

    fin.close();
    std::cout << "Read from file " << fileName << std::endl;
    std::cout << "  DIM=" << DIM << " n=" << n << " steps=" << steps << " dt=" << dt << std::endl;
}

inline void initAndSaveRandom() {
    int n = 100;
    int steps = 10;
    double dt = 0.01;

    std::vector<double> masses(n);
    std::vector<Vec> positions(n);
    std::vector<Vec> velocities(n);

    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> posd(-1e2, 1e2);
    std::uniform_real_distribution<double> massd(1e3, 1e3);
    std::uniform_real_distribution<double> veld(-1e1, 1e1);
    for (int i = 0; i < n; ++i) {
        masses[i] = massd(rng);
        for (int j = 0; j < DIM; j++) {
            positions[i][j] = posd(rng);
            velocities[i][j] = veld(rng);
        }
    }

    // save to file
    saveToFile("test1.in.out", n, steps, dt, masses, positions, velocities);
}

inline void runSerial() {

    // BELOW CODE: for generating random test1.in.out
    // initAndSaveRandom();
    // return 0;

    int n;
    int steps;
    double dt;

    std::vector<double> masses;
    std::vector<Vec> positions;
    std::vector<Vec> velocities;
    std::vector<Vec> forces;

    readFromFile("test1.in.out", n, steps, dt, masses, positions, velocities);

    forces.resize(n);
    for (int i = 0; i < n; i++) for (int j = 0; j < DIM; j++)
        forces[i][j] = 0.;

    for (int step = 0; step < steps; step++) {
        // make forces equal to zero
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < DIM; j++) {
                forces[i][j] = 0.;
            }
        }

        // compute all forces
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < DIM; j++) {
                forces[i][j] = 0.;
            }

            for (int k = 0; k < n; k++) {
                if (i == k) continue;
                Vec force = force_on_p1(positions[i], positions[k], masses[i], masses[k]);
                for (int j = 0; j < DIM; j++) {
                    forces[i][j] += force[j];
                }
            }
        }

        // apply forces, i.e. calc new positions and velocities
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < DIM; j++) {
                positions[i][j] += dt * velocities[i][j];
                velocities[i][j] += dt / masses[i] * forces[i][j];
            }
        }
    
        // save for testing
        saveToFile("test1." + std::to_string(step) + ".out", n, steps, dt,
            masses, positions, velocities);
    }

}

int main(int argc, char** argv) {

}