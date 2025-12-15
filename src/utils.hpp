#pragma once

#include <assert.h>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "body.hpp"

namespace utils {

template <int DIM>
void saveToStream(std::ostream &out, int steps, double dt,
                  const bodies<DIM>& bodies, bool saveVelocities = true) {
  out << DIM << "\n";
  out << bodies.globalSize() << " " << steps << " " << dt << "\n";
  for (size_t i = 0; i < bodies.globalSize(); i++) {
    out << bodies.global(i).mass() << " ";
    for (int j = 0; j < DIM; j++)
      out << bodies.global(i).pos()[j] << " ";
    if (saveVelocities) {
      for (int j = 0; j < DIM; j++)
        out << bodies.global(i).vel()[j] << " ";
    }
    out << "\n";
  }
}

template <int DIM>
void readFromStream(std::istream &in, int &steps, double &dt,
                    bodies<DIM>& bodies, bool readVelocities = true) {
  int dimInFile;
  in >> dimInFile;
  assert(dimInFile == DIM);

  int n;

  in >> n >> steps >> dt;

  bodies.resize(n, n, 0);

  for (int i = 0; i < n; i++) {
    in >> bodies.global(i).mass();
    for (int j = 0; j < DIM; j++)
      in >> bodies.global(i).pos()[j];
    if (readVelocities) {
      for (int j = 0; j < DIM; j++)
        in >> bodies.global(i).vel()[j];
    }
  }
}

template <int DIM>
void saveToFile(const std::string &fileName, int steps, double dt,
                const bodies<DIM>& bodies, bool saveVelocities = true) {
  std::ofstream fout(fileName);

  saveToStream(fout, steps, dt, bodies, saveVelocities);

  fout.flush();
  fout.close();

  std::cout << "Saved to file " << fileName << std::endl;
}

template <int DIM>
void readFromFile(const std::string &fileName, int &steps, double &dt,
                  bodies<DIM>& bodies, bool readVelocities = true) {
  std::ifstream fin(fileName);

  readFromStream(fin, steps, dt, bodies, readVelocities);

  fin.close();
  std::cout << "Read from file " << fileName << std::endl;
  std::cout << "  DIM=" << DIM << " n=" << bodies.globalSize() << " steps=" << steps
            << " dt=" << dt << std::endl;
}

template <int DIM>
void generateRandomToFile(const std::string &filename, int n = 100,
                          int steps = 10, double dt = 0.01, int seed = 42) {
  constexpr int MAX_POS_COMPONENT = 100;
  constexpr int MAX_VEL_COMPONENT = 10;
  constexpr int MIN_MASS = 1;
  constexpr int MAX_MASS = 100;

  std::vector<double> masses(n);
  std::vector<Vec<DIM>> positions(n);
  std::vector<Vec<DIM>> velocities(n);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> posd(-MAX_POS_COMPONENT,
                                              MAX_POS_COMPONENT);
  std::uniform_real_distribution<double> massd(MIN_MASS, MAX_MASS);
  std::uniform_real_distribution<double> veld(-MAX_VEL_COMPONENT,
                                              MAX_VEL_COMPONENT);
  for (int i = 0; i < n; ++i) {
    masses[i] = massd(rng);
    for (int j = 0; j < DIM; j++) {
      positions[i][j] = posd(rng);
      velocities[i][j] = veld(rng);
    }
  }

  saveToFile(filename, n, steps, dt, masses, positions, velocities);
}

} // namespace utils
