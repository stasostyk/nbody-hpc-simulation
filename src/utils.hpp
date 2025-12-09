#pragma once

#include <iostream>
#include <vector>

#include "core.hpp"

namespace utils {

void saveToStream(std::ostream &out, int n, int steps, double dt,
                  const std::vector<double> &m, const std::vector<Vec> &p,
                  const std::vector<Vec> &v, bool saveVelocities = true);

void readFromStream(std::istream &in, int &n, int &steps, double &dt,
                    std::vector<double> &m, std::vector<Vec> &p,
                    std::vector<Vec> &v, bool readVelocities = true);

void saveToFile(const std::string &fileName, int n, int steps, double dt,
                const std::vector<double> &m, const std::vector<Vec> &p,
                const std::vector<Vec> &v, bool saveVelocities = true);

void readFromFile(const std::string &fileName, int &n, int &steps, double &dt,
                  std::vector<double> &m, std::vector<Vec> &p,
                  std::vector<Vec> &v, bool readVelocities = true);

void generateRandomToFile(const std::string &filename, int n = 100,
                          int steps = 10, double dt = 0.01, int seed = 42);

} // namespace utils
