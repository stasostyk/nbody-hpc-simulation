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
                  const Bodies<DIM, EmptyAttributes>& bodies, bool saveVelocities = true) {
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
                    Bodies<DIM, EmptyAttributes>& bodies, bool readVelocities = true) {
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
                const Bodies<DIM, EmptyAttributes>& bodies, bool saveVelocities = true) {
  std::ofstream fout(fileName);

  saveToStream(fout, steps, dt, bodies, saveVelocities);

  fout.flush();
  fout.close();

  std::cout << "Saved to file " << fileName << std::endl;
}

template <int DIM>
void readFromFile(const std::string &fileName, int &steps, double &dt,
                  Bodies<DIM, EmptyAttributes>& bodies, bool readVelocities = true) {
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    throw std::runtime_error("File not found: " + fileName + ". Ensure that it exists then try again.");
  }

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

  Bodies<DIM, EmptyAttributes> bodies;
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

  saveToFile(filename, steps, dt, bodies);
}

template<int DIM>
void compareOutputs(const std::string& prefixA = "test1", const std::string& prefixB = "test1-MPI", int steps = 10) {
    constexpr double tol = 1e-10;

    for (int step = 0; step < steps; step++) {
        int stepsA, stepsB;
        double dtA, dtB;
        Bodies<DIM, EmptyAttributes> bodiesA, bodiesB;

        std::string filenameSuffix = + "." + std::to_string(step) + ".out";
        std::string filenameA = prefixA + filenameSuffix;
        std::string filenameB = prefixB + filenameSuffix;

        utils::readFromFile(filenameA, stepsA, dtA, bodiesA);
        utils::readFromFile(filenameB, stepsB, dtB, bodiesB);

        assert(steps == stepsA && stepsA == stepsB);
        assert(bodiesA.localSize() == bodiesB.localSize());
        assert(fabs(dtA - dtB) < tol);
        
        for (size_t i = 0; i < bodiesA.localSize(); i++) {
            assert(fabs(bodiesA.local(i).mass() == bodiesB.local(i).mass()) < tol);
            for (int j = 0; j < DIM; j++) {
                assert(fabs(bodiesA.local(i).pos()[j] == bodiesB.local(i).pos()[j]) < tol);
            }
        }
    }

    std::cout << "CHECK FINISHED, EVERYTHING IS FINE" << std::endl;
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


template<int DIM>
void compareOutputsSingleFile(std::string filename1, std::string filename2) {
  std::cout << "===============================================" << std::endl;
  std::cout << "Comparing files " << filename1 << " and " << filename2 << std::endl;

  double e = 1e-6;

  int stepsA, stepsB;
  double dtA, dtB;
  Bodies<DIM, EmptyAttributes> bodiesA, bodiesB;

  utils::readFromFile<DIM>(filename1, stepsA, dtA, bodiesA);
  utils::readFromFile<DIM>(filename2, stepsB, dtB, bodiesB);
  int n1 = bodiesA.globalSize();
  int n2 = bodiesB.globalSize();

  check(stepsA == stepsB, "steps");
  check(n1 == n2, "n");
  check(fabs(dtA - dtB) < e, "dt");
  
  for (int i = 0; i < n1; i++) {
      // check(fabs(masses1[i]- masses2[i]) < e, "masses");
      checkDoubles(bodiesA.global(i).mass(), bodiesB.global(i).mass(), "masses " + std::to_string(i));
      for (int j = 0; j < DIM; j++) {
          // check(fabs(positions1[i][j] - positions2[i][j]) < e, "positions");
          checkDoubles(bodiesA.global(i).pos()[j], bodiesB.global(i).pos()[j], "positions " + std::to_string(i));
      }
  }

  std::cout << "CHECK FINISHED, EVERYTHING IS FINE" << std::endl;
  std::cout << "===============================================" << std::endl;

}


} // namespace utils
