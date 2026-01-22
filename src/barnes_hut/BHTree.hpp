#pragma once

#include <array>
#include <vector>
#include <memory>
#include <iostream>

#include "bodies.hpp"
#include "../body.hpp"
#include "../vec.hpp"
#include "../forces/func.hpp"

const double boundCheckEps = 1e-10;
const double distanceEps = 1e-10;

template<int dim>
class Box {
public:
    Box(const Box<dim> &other) : corner1(other.corner1), corner2(other.corner2) {}

    Box(const Vec<dim> &_corner1, const Vec<dim> &_corner2) 
    : corner1(_corner1), corner2(_corner2) {} 

    bool isPointInside(const Vec<dim> &position) const {
        for (int dimId = 0; dimId < dim; dimId++) {
            const double val = position[dimId];
            if (val < corner1[dimId] - boundCheckEps) {
                return false;
            }
            if (val >= corner2[dimId]) {
                return false;
            }
        }

        return true;
    }

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
    }

    ~BHTree() = default;
    
    void printInfo() const {
        std::cout << "=========================================" << std::endl;
        std::cout << "BH Tree. dim = " << dim << std::endl;
        root->printInfoRecursively();
        std::cout << "=========================================" << std::endl;
    }

    Vec<dim> calculateForce(const Body<dim>& body, double theta, const forces::force<dim, Attributes> &force) const {
        return root->calculateForce(body, theta, force);
    }

private:


    class Node {
    public:
        Node (const Box<dim> &_bounds, const int &_depth)
        : bounds(_bounds), depth(_depth), isLeaf(true), isEmptyNode(true)
        {
            centerOfMass.mass = 0.;
            centerOfMass.position = 0.;
            centerOfMass.velocity = 0.;

            // Using one side length as length l used for Barnes Hut threshold s/d<theta.
            const double length = bounds.corner2[0] - bounds.corner1[0];
            lengthSqr = length * length;

            activeChildren.reserve(childrenCnt);
        }

        void insertBody(const Body<dim> &body) {
            if (isEmptyNode) {
                isEmptyNode = false;
                centerOfMass.combineWith(body);
                return;
            }

            centerOfMass.combineWith(body);
            insertToChild(body);
        }

        void printInfoRecursively() const {
            std::cout << "---" << std::endl;
            std::cout << "Node, depth = " << depth << std::endl;
            std::cout << " bounds from: ";
            bounds.corner1.print();
            std::cout << " to: ";
            bounds.corner2.print();
            std::cout << std::endl;
            std::cout << " mass: " << centerOfMass.mass << ", centerOfMass: ";
            centerOfMass.position.print();
            std::cout << std::endl;

            if (isLeaf) return;
            for (int childId = 0; childId < childrenCnt; childId++) {
                children[childId]->printInfoRecursively();
            }
        }
        
        Vec<dim> calculateForce(const Body<dim>& body, double theta, const forces::force<dim, Attributes> &force) const {
            if (isLeaf) {
                // in the leaf node there should be max one body, so center of mass is the body
                if (!isEmptyNode && centerOfMass.bodyId != body.bodyId) {
                    return force(
                        bodyCopy<dim, EmptyAttributes>(body.position, Vec<dim> (0.), body.mass, EmptyAttributes()),
                        bodyCopy<dim, EmptyAttributes>(centerOfMass.position, Vec<dim> (0.), centerOfMass.mass, EmptyAttributes())
                    );
                }

                return Vec<dim>(0.);
            }

            if (isEmptyNode) {
                return Vec<dim>(0.);
            }

            const double distanceSqr = (centerOfMass.position - body.position).normSquared();

            const double bhVal = lengthSqr / (distanceSqr + distanceEps*distanceEps);

            // Barnes-Hut criterion
            if (bhVal < theta) {
                // Approximating far away bodies as one with computed center of mass.
                return force(
                    bodyCopy<dim, EmptyAttributes>(body.position, Vec<dim> (0.), body.mass, EmptyAttributes()),
                    bodyCopy<dim, EmptyAttributes>(centerOfMass.position, Vec<dim> (0.), centerOfMass.mass, EmptyAttributes())
                );
            }


            Vec<dim> totalForce(0.);
            for (int childId : activeChildren) {
                totalForce += children[childId]->calculateForce(body, theta, force);
            }

            return totalForce;
        }


    private:
        static constexpr int childrenCnt = 1 << dim;

        Box<dim> bounds;
        const int depth;
        double lengthSqr;
        bool isLeaf;
        bool isEmptyNode; // doesn't have any bodies inside
        std::array<std::unique_ptr<Node>, childrenCnt> children;
        std::vector<int> activeChildren;

        // Center of mass of all bodies in the node.
        // Has averaged position and velocity.
        Body<dim> centerOfMass;

        void insertToChild(const Body<dim> &body) {
            int childId = getChildIndex(body.position);

            if (!children[childId]) {
                initChild(childId);
            }

            children[childId]->insertBody(body);
        }


        int getChildIndex(const Vec<dim> &position) const {
            const Vec<dim> midPoint = (bounds.corner1 + bounds.corner2) * 0.5;
            int index = 0;
            for (int dimId = 0; dimId < dim; dimId++) {
                if (position[dimId] >= midPoint[dimId] - boundCheckEps) {
                    index |= (1 << dimId);
                }
            }
            return index;
        }
        
        void initChild(int childId) {
            Vec<dim> corner1;
            Vec<dim> corner2;

            const Vec<dim> midPoint = (bounds.corner1 + bounds.corner2) * 0.5;

            for (int dimId = 0; dimId < dim; dimId++) {
                if ((childId & (1 << dimId)) == 0) {
                    corner1[dimId] = bounds.corner1[dimId];
                    corner2[dimId] = midPoint[dimId];
                } else {
                    corner1[dimId] = midPoint[dimId];
                    corner2[dimId] = bounds.corner2[dimId];
                }
            }

            const Box<dim> childBounds(corner1, corner2);
            children[childId] = std::make_unique<Node>(childBounds, depth + 1);
            activeChildren.push_back(childId);  
        }


    };

    std::unique_ptr<Node> root;
    const Box<dim> universeBounds;
};
