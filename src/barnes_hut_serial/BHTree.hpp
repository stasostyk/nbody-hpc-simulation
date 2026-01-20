#pragma once

#include <array>
#include <vector>
#include <memory>
#include <iostream>

#include "../bodies.hpp"
#include "../vec.hpp"

const double boundCheckEps = 1e-10;

template<int dim>
class Box {
public:
    Box(const Vec<dim> &_corner1, const Vec<dim> &_corner2) 
    : corner1(_corner1), corner2(_corner2) {} 

    bool isPointInside(const Vec<dim> &position) const;

    Vec<dim> corner1;
    Vec<dim> corner2;
};

// Implementation of Barnes-Hut tree.
// In 2D works as Quadtree, in 3D as Octree.
template <int dim>
class BHTree {
public:

    BHTree(const Box<dim> &_universeBounds, const double _theta, const std::vector<Body<dim>> _bodies)
    : universeBounds(_universeBounds), theta(_theta) {
        root = std::make_unique<Node>(_universeBounds, 0);
        for (const Body<dim> &body : _bodies) {
            root->insertBody(body);
        }
        root->calculateCentersOfMassRecursively();
    }

    ~BHTree() = default;
    void printInfo() const;
    Vec<dim> calculateForce(const Body<dim>& b) const;

private:

    class Node {
    public:
        Node (const Box<dim> &_bounds, const int &_depth)
        : bounds(_bounds), depth(_depth), isLeaf(true),
          totalMass(0.), centerOfMass(0.)
        {}

        void insertBody(const Body<dim> &body);
        void printInfoRecursively() const;
        void calculateCentersOfMassRecursively();
        Vec<dim> calculateForce(const Body<dim> &body) const;

    private:
        static constexpr int childrenCnt = 1 << dim;

        void splitNode();
        void insertToChild(const Body<dim> &body);
        void calculateCenterOfMass();

        Box<dim> bounds;
        const int depth;
        bool isLeaf;
        std::array<std::unique_ptr<Node>, childrenCnt> children;
        std::vector<Body<dim>> bodies; // TODO probably no need to store all bodies

        // Center of mass of all bodies in the node
        double totalMass;
        Vec<dim> centerOfMass;

    };

    std::unique_ptr<Node> root;
    const Box<dim> universeBounds;
    const double theta; // BH threshold for s/d
};
