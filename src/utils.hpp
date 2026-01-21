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
  if (!fin.is_open()) {
    throw std::runtime_error("File not found: " + fileName + ". Ensure that it exists then try again.");
  }

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

void check(bool cond, std::string msg) {
  if (cond) return;
  std::cout << "CONDITION FAILED, msg: " << msg << std::endl;
  exit(1);
}

void checkDoubles(double d1, double d2, std::string desc, double harshEps = 1e-2, double softEps = 1e-4) {
  double diff = fabs(d1 - d2);
  if (diff > softEps) {
    std::cout << "Error large for " << desc << ", val1=" << d1 << ", val2=" << d2 << std::endl;
  }
  check(diff < harshEps, desc);
}

template <int DIM>
void compareOutputs(std::string filename1, std::string filename2) {
  std::cout << "===============================================" << std::endl;
  std::cout << "Comparing files " << filename1 << " and " << filename2 << std::endl;

  double e = 1e-6;

  int n1, n2;
  int steps1, steps2;
  double dt1, dt2;

  std::vector<double> masses1, masses2;
  std::vector<Vec<DIM>> positions1, positions2;

  utils::readFromFile(filename1, n1, steps1, dt1, 
                masses1, positions1, positions1, false);
  utils::readFromFile(filename2, n2, steps2, dt2, 
                masses2, positions2, positions2, false);

  check(steps1 == steps2, "steps");
  check(n1 == n2, "n");
  check(fabs(dt1 - dt2) < e, "dt");
  
  for (int i = 0; i < n1; i++) {
      // check(fabs(masses1[i]- masses2[i]) < e, "masses");
      checkDoubles(masses1[i], masses2[i], "masses " + std::to_string(i));
      for (int j = 0; j < DIM; j++) {
          // check(fabs(positions1[i][j] - positions2[i][j]) < e, "positions");
          checkDoubles(positions1[i][j], positions2[i][j], "positions " + std::to_string(i));
      }
  }

  std::cout << "CHECK FINISHED, EVERYTHING IS FINE" << std::endl;
  std::cout << "===============================================" << std::endl;

}


} // namespace utils
