#include "utils.hpp"
#include <cstdlib>
#include <iostream>

constexpr int dim = 2;

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " [filename1] [filename2]\n";
        return 1;
    }

    std::string filename1 = argv[1];
    std::string filename2 = argv[2];

    utils::compareOutputsSingleFile<dim>(filename1, filename2);

    return 0;
}