#pragma once

#include <assert.h>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "vec.hpp"

namespace utils {

template <int DIM>
void saveToStream(std::ostream &out, int n, int steps, double dt,
                  const std::vector<double> &m, const std::vector<Vec<DIM>> &p,
                  const std::vector<Vec<DIM>> &v, bool saveVelocities = true) {
  out << DIM << "\n";
  out << n << " " << steps << " " << dt << "\n";
  for (int i = 0; i < n; i++) {
    out << m[i] << " ";
    for (int j = 0; j < DIM; j++)
      out << p[i][j] << " ";
    if (saveVelocities) {
      for (int j = 0; j < DIM; j++)
        out << v[i][j] << " ";
    }
    out << "\n";
  }
}

template <int DIM>
void readFromStream(std::istream &in, int &n, int &steps, double &dt,
                    std::vector<double> &m, std::vector<Vec<DIM>> &p,
                    std::vector<Vec<DIM>> &v, bool readVelocities = true) {
  int dimInFile;
  in >> dimInFile;
  assert(dimInFile == DIM);

  in >> n >> steps >> dt;

  m.resize(n);
  p.resize(n);
  v.resize(n);

  for (int i = 0; i < n; i++) {
    in >> m[i];
    for (int j = 0; j < DIM; j++)
      in >> p[i][j];
    if (readVelocities) {
      for (int j = 0; j < DIM; j++)
        in >> v[i][j];
    }
  }
}

template <int DIM>
void saveToFile(const std::string &fileName, int n, int steps, double dt,
                const std::vector<double> &m, const std::vector<Vec<DIM>> &p,
                const std::vector<Vec<DIM>> &v, bool saveVelocities = true) {
  std::ofstream fout(fileName);

  saveToStream(fout, n, steps, dt, m, p, v, saveVelocities);

  fout.flush();
  fout.close();

  std::cout << "Saved to file " << fileName << std::endl;
}

template <int DIM>
void readFromFile(const std::string &fileName, int &n, int &steps, double &dt,
                  std::vector<double> &m, std::vector<Vec<DIM>> &p,
                  std::vector<Vec<DIM>> &v, bool readVelocities = true) {
  std::ifstream fin(fileName);

  readFromStream(fin, n, steps, dt, m, p, v, readVelocities);

  fin.close();
  std::cout << "Read from file " << fileName << std::endl;
  std::cout << "  DIM=" << DIM << " n=" << n << " steps=" << steps
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
