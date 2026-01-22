#pragma once

#include "vec.hpp"

template <int DIM, typename Attributes> 
struct body;
template <int DIM, typename Attributes> 
struct bodyView;

struct EmptyAttributes {};

// Store bodies as structure of arrays
// View bodies as array of structure
template <int DIM, typename Attributes = EmptyAttributes> 
struct Bodies {
private:
  size_t _globalSize = 0;
  size_t _localSize = 0;
  size_t _localOffset = 0;
  std::vector<bodyView<DIM, Attributes>> _views;

public:
  std::vector<Vec<DIM>> position;
  std::vector<Vec<DIM>> velocity;
  std::vector<double> mass;
  std::vector<Attributes> attributes;

  void resize(int globalSize, int localSize, int localOffset) {
    _globalSize = globalSize;
    _localSize = localSize;
    _localOffset = localOffset;

    position.resize(globalSize);
    velocity.resize(globalSize);
    mass.resize(globalSize);
    attributes.resize(globalSize);

    _views.clear();
    for (int i = 0; i < globalSize; i++) {
      _views.push_back(bodyView<DIM, Attributes>(*this, i));
    }
  }
  body<DIM, Attributes> &global(int index) { return _views[index]; }
  const body<DIM, Attributes> &global(int index) const { return _views[index]; }
  size_t globalSize() const { return _globalSize; }
  body<DIM, Attributes> &local(int index) { return _views[_localOffset + index]; }
  const body<DIM, Attributes> &local(int index) const {
    return _views[_localOffset + index];
  }
  size_t localSize() const { return _localSize; }
  size_t localOffset() const { return _localOffset; }
};

template <int DIM, typename Attributes> 
struct body {
public:
  virtual Vec<DIM> &pos() = 0;
  virtual const Vec<DIM> &pos() const = 0;
  virtual Vec<DIM> &vel() = 0;
  virtual const Vec<DIM> &vel() const = 0;
  virtual double &mass() = 0;
  virtual const double &mass() const = 0;
  virtual Attributes &attributes() = 0;
  virtual const Attributes &attributes() const = 0;
};

template <int DIM, typename Attributes> 
struct bodyView : public body<DIM, Attributes> {
private:
  Bodies<DIM, Attributes> &_bodies;

public:
  bodyView(Bodies<DIM, Attributes> &bodies, int index) : _bodies(bodies), index(index) {}
  const int index;
  Vec<DIM> &pos() override { return _bodies.position[index]; }
  const Vec<DIM> &pos() const override { return _bodies.position[index]; }
  Vec<DIM> &vel() override { return _bodies.velocity[index]; }
  const Vec<DIM> &vel() const override { return _bodies.velocity[index]; }
  double &mass() override { return _bodies.mass[index]; }
  const double &mass() const override { return _bodies.mass[index]; }
  Attributes &attributes() { return _bodies.attributes[index]; }
  const Attributes &attributes() const { return _bodies.attributes[index]; }
};

template <int DIM, typename Attributes> 
struct bodyCopy : public body<DIM, Attributes> {
private:
  Vec<DIM> _pos;
  Vec<DIM> _vel;
  double _mass;
  Attributes _attributes;

public:
  bodyCopy(Vec<DIM> pos, Vec<DIM> vel, double mass, Attributes attributes) {
      std::copy_n(pos.begin(), DIM, _pos.begin());
      std::copy_n(vel.begin(), DIM, _vel.begin());
      _mass = mass;
      _attributes = attributes;
  }
  Vec<DIM> &pos() override { return _pos; }
  const Vec<DIM> &pos() const override { return _pos; }
  Vec<DIM> &vel() override { return _vel; }
  const Vec<DIM> &vel() const override { return _vel; }
  double &mass() override { return _mass; }
  const double &mass() const override { return _mass; }
  Attributes &attributes() { return _attributes; }
  const Attributes &attributes() const { return _attributes; }
};
