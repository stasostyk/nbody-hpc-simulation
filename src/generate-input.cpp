#include "utils.hpp"
#include <cstdlib>
#include <iostream>

constexpr int DIM = 3;

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " [n] [steps] [dt] [seed]\n";
        return 1;
    }

    std::string filename = "test1.in.out";
    int n     = (argc > 1) ? std::atoi(argv[1]) : 100;
    int steps = (argc > 2) ? std::atoi(argv[2]) : 10;
    double dt = (argc > 3) ? std::atof(argv[3]) : 0.01;
    int seed  = (argc > 4) ? std::atoi(argv[4]) : 42;

    utils::generateRandomToFile<DIM>(filename, n, steps, dt, seed);

    return 0;
}