#pragma once

#include <assert.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "body.hpp"

namespace utils {

template <int DIM>
void saveToStream(std::ostream &out, int n, double total_time, double dt_max,
                  const bodies<DIM, EmptyAttributes> &bodies,
                  bool saveVelocities = true) {
  out << DIM << "\n";
  out << n << " " << total_time << " " << dt_max << "\n";
  for (int i = 0; i < n; i++) {
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
void readFromStream(std::istream &in, int &n, double &total_time,
                    double &dt_max, bodies<DIM, EmptyAttributes> &bodies,
                    bool readVelocities = true) {
  int dimInFile;
  in >> dimInFile;
  assert(dimInFile == DIM);

  in >> n >> total_time >> dt_max;

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
void saveToFile(const std::string &fileName, int n, double total_time,
                double dt_max, const bodies<DIM, EmptyAttributes> &bodies,
                bool saveVelocities = true) {
  std::ofstream fout(fileName);

  saveToStream(fout, n, total_time, dt_max, bodies, saveVelocities);

  fout.flush();
  fout.close();

  std::cout << "Saved to file " << fileName << std::endl;
}

template <int DIM>
void readFromFile(const std::string &fileName, int &n, double &total_time,
                  double &dt_max, bodies<DIM, EmptyAttributes> &bodies,
                  bool readVelocities = true) {
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    throw std::runtime_error("File not found: " + fileName +
                             ". Ensure that it exists then try again.");
  }

  readFromStream(fin, n, total_time, dt_max, bodies, readVelocities);

  fin.close();
  std::cout << "Read from file " << fileName << std::endl;
  std::cout << "  DIM=" << DIM << " n=" << n << " T_end=" << total_time
            << " dt_max=" << dt_max << std::endl;
}

template <int DIM>
void generateRandomToFile(const std::string &filename, int n = 100,
                          int steps = 10, double dt = 0.01, int seed = 42) {
  constexpr int MAX_POS_COMPONENT = 100;
  constexpr int MAX_VEL_COMPONENT = 10;
  constexpr int MIN_MASS = 1;
  constexpr int MAX_MASS = 100;

  bodies<DIM, EmptyAttributes> bodies;
  bodies.resize(n, n, 0);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> posd(-MAX_POS_COMPONENT,
                                              MAX_POS_COMPONENT);
  std::uniform_real_distribution<double> massd(MIN_MASS, MAX_MASS);
  std::uniform_real_distribution<double> veld(-MAX_VEL_COMPONENT,
                                              MAX_VEL_COMPONENT);
  for (int i = 0; i < n; ++i) {
    bodies.mass[i] = massd(rng);
    for (int j = 0; j < DIM; j++) {
      bodies.position[i][j] = posd(rng);
      bodies.velocity[i][j] = veld(rng);
    }
  }

  double total_time = static_cast<double>(steps) * dt;
  saveToFile(filename, n, total_time, dt, bodies);
}

template <int DIM>
void compareOutputs(const std::string &prefixA = "test1",
                    const std::string &prefixB = "test1-MPI", int steps = 10) {
  constexpr double tol = 1e-10;

  for (int step = 0; step < steps; step++) {
    int nA, nB;
    double total_timeA, total_timeB;
    double dt_maxA, dt_maxB;
    bodies<DIM, EmptyAttributes> bodiesA, bodiesB;

    std::string filenameSuffix = "." + std::to_string(step) + ".out";
    std::string filenameA = prefixA + filenameSuffix;
    std::string filenameB = prefixB + filenameSuffix;

    utils::readFromFile(filenameA, nA, total_timeA, dt_maxA, bodiesA);
    utils::readFromFile(filenameB, nB, total_timeB, dt_maxB, bodiesB);

    assert(nA == nB);
    assert(fabs(total_timeA - total_timeB) < tol);
    assert(fabs(dt_maxA - dt_maxB) < tol);

    for (int i = 0; i < nA; i++) {
      assert(fabs(bodiesA.global(i).mass() - bodiesB.global(i).mass()) < tol);
      for (int j = 0; j < DIM; j++) {
        assert(fabs(bodiesA.global(i).pos()[j] - bodiesB.global(i).pos()[j]) <
               tol);
      }
    }
  }

  std::cout << "CHECK FINISHED, EVERYTHING IS FINE" << std::endl;
}

} // namespace utils
