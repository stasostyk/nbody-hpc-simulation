#include "utils.hpp"
#include <cstdlib>
#include <iostream>

constexpr int DIM = 3;

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " [filename] [n] [steps] [dt] [seed]\n";
        return 1;
    }

    std::string filename = argv[1];
    int n     = (argc > 2) ? std::atoi(argv[2]) : 100;
    int steps = (argc > 3) ? std::atoi(argv[3]) : 10;
    double dt = (argc > 4) ? std::atof(argv[4]) : 0.01;
    int seed  = (argc > 5) ? std::atoi(argv[5]) : 42;

    utils::generateRandomToFile<DIM>(filename, n, steps, dt, seed);

    return 0;
}