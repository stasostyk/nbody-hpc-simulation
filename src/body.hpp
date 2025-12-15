#pragma once

#include "vec.hpp"

template <int DIM> struct body;

// Store bodies as structure of arrays
// View bodies as array of structure
template <int DIM> struct bodies {
private:
  size_t _globalSize = 0;
  size_t _localSize = 0;
  size_t _localOffset = 0;
  std::vector<body<DIM>> _views;

public:
  std::vector<Vec<DIM>> position;
  std::vector<Vec<DIM>> velocity;
  std::vector<double> mass;

  void resize(int globalSize, int localSize, int localOffset) {
    _globalSize = globalSize;
    _localSize = localSize;
    _localOffset = localOffset;

    position.resize(globalSize);
    velocity.resize(globalSize);
    mass.resize(globalSize);

    _views.clear();
    for (int i = 0; i < globalSize; i++) {
      _views.push_back(body<DIM>(*this, i));
    }
  }
  body<DIM> &global(int index) { return _views[index]; }
  const body<DIM> &global(int index) const { return _views[index]; }
  size_t globalSize() const { return _globalSize; }
  body<DIM> &local(int index) { return _views[_localOffset + index]; }
  const body<DIM> &local(int index) const {
    return _views[_localOffset + index];
  }
  size_t localSize() const { return _localSize; }
  size_t localOffset() const { return _localOffset; }
};

template <int DIM> struct body {
private:
  bodies<DIM>& _bodies;

public:
  body(bodies<DIM> &bodies, int index) : _bodies(bodies), index(index) {}
  const int index;
  Vec<DIM> &pos() { return _bodies.position[index]; }
  const Vec<DIM> &pos() const { return _bodies.position[index]; }
  Vec<DIM> &vel() { return _bodies.velocity[index]; }
  const Vec<DIM> &vel() const { return _bodies.velocity[index]; }
  double &mass() { return _bodies.mass[index]; }
  const double &mass() const { return _bodies.mass[index]; }
};
