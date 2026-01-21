#pragma once

#include <array>
#include <vector>
#include <memory>
#include <iostream>

#include "../bodies.hpp"
#include "../body.hpp"
#include "../vec.hpp"
#include "../forces/func.hpp"

const double boundCheckEps = 1e-10;
const double distanceEps = 1e-10;

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
template <int dim, typename Attributes>
class BHTree {
public:

    BHTree(
        const Box<dim> &_universeBounds, 
        const std::vector<Body<dim>> &_bodies
    ) : universeBounds(_universeBounds) {
        root = std::make_unique<Node>(_universeBounds, 0);
        for (const Body<dim> &body : _bodies) {
            root->insertBody(body);
        }
        // root->calculateCentersOfMassRecursively();
        // std::cout << "MAX DEPTH: " << root->getMaxDepth() << std::endl;
    }

    ~BHTree() = default;
    void printInfo() const;
    Vec<dim> calculateForce(const Body<dim>& body, double theta, const forces::force<dim, Attributes> &force) const;

private:

    class Node {
    public:
        Node (const Box<dim> &_bounds, const int &_depth)
        : bounds(_bounds), depth(_depth), isLeaf(true), isEmptyNode(true)
        {
            centerOfMass.mass = 0.;
            centerOfMass.position = 0.;
            centerOfMass.velocity = 0.;

            // TODO now assumes the bounds are square/cube
            const double length = bounds.corner2[0] - bounds.corner1[0];
            lengthSqr = length * length;

            activeChildren.reserve(childrenCnt);
        }

        void insertBody(const Body<dim> &body);
        void printInfoRecursively() const;
        // void calculateCentersOfMassRecursively();
        Vec<dim> calculateForce(const Body<dim> &body, double theta, const forces::force<dim, Attributes> &force) const;
        int getMaxDepth();

    private:
        static constexpr int childrenCnt = 1 << dim;

        void splitNode();
        void insertToChild(const Body<dim> &body);
        int getChildIndex(const Vec<dim> &position) const;
        void initChild(int childId);
        // void calculateCenterOfMass();
        
        Box<dim> bounds;
        const int depth;
        double lengthSqr;
        bool isLeaf;
        bool isEmptyNode; // doesn't have any bodies inside
        std::array<std::unique_ptr<Node>, childrenCnt> children;
        std::vector<int> activeChildren;
        // std::vector<Body<dim>> bodies; // TODO probably no need to store all bodies

        // Center of mass of all bodies in the node.
        // Has averaged position and velocity.
        Body<dim> centerOfMass;

    };

    std::unique_ptr<Node> root;
    const Box<dim> universeBounds;
};
