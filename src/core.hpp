#include <array>

// Dimensions of the problem (2D or 3D)
constexpr int DIM = 2;

// Gravitational constant
constexpr double G = 6.673e-11;

// Vector in our DIM-dimensional space
using Vec = std::array<double, DIM>;
