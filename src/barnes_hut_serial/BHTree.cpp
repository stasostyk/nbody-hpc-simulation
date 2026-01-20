#include "BHTree.hpp"

template<int dim>
void BHTree<dim>::printInfo() const {
    std::cout << "=========================================" << std::endl;
    std::cout << "BH Tree. dim = " << dim << std::endl;
    root->printInfoRecursively();
    std::cout << "=========================================" << std::endl;
}

template<int dim>
bool Box<dim>::isPointInside(const Vec<dim> &position) const {
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

template<int dim>
void BHTree<dim>::Node::insertBody(const Body<dim> &body) {
    if (bodies.empty()) {
        bodies.push_back(body);
        return;
    }

    if (isLeaf) {
        splitNode(); // Make sure there are children.
    }

    bodies.push_back(body);
    insertToChild(body);
}

template<int dim>
void BHTree<dim>::Node::insertToChild(const Body<dim> &body) {
    // Find to which child node the new body belongs
    for (int childId = 0; childId < childrenCnt; childId++) {
        Vec<dim> bodyPos(body.position);
        if (children[childId]->bounds.isPointInside(bodyPos)) {
            children[childId]->insertBody(body);
            return;
        }
    }
}

template<int dim>
void BHTree<dim>::Node::splitNode() {
    if (!isLeaf) {
        // Already split.
        return;
    }
    isLeaf = false;

    for (int bitmask = 0; bitmask < (int)childrenCnt; bitmask++) {
        Vec<dim> corner1;
        Vec<dim> corner2;

        const Vec<dim> midPoint = (bounds.corner1 + bounds.corner2) / 2.;

        for (int dimId = 0; dimId < dim; dimId++) {
            if ((bitmask & (1 << dimId)) == 0) {
                // bitmask[dimId] is turned off
                corner1[dimId] = bounds.corner1[dimId];
                corner2[dimId] = midPoint[dimId];
            } else {
                // bitmask[dimId] is turned on
                corner1[dimId] = midPoint[dimId];
                corner2[dimId] = bounds.corner2[dimId];
            }
        }

        const Box<dim> childBounds(corner1, corner2);
        children[bitmask] = std::make_unique<Node>(childBounds, depth + 1);
    }

    for (const Body<dim> &body : bodies) {
        insertToChild(body);
    }
}


template<int dim>
void BHTree<dim>::Node::printInfoRecursively() const {
    std::cout << "---" << std::endl;
    std::cout << "Node, depth = " << depth << std::endl;
    std::cout << " bounds from: ";
    bounds.corner1.print();
    std::cout << " to: ";
    bounds.corner2.print();
    std::cout << std::endl;
    std::cout << " mass: " << totalMass << ", centerOfMass: ";
    centerOfMass.print();
    std::cout << std::endl;

    std::cout << " bodies (positions): " << std::endl;
    for (const Body<dim> &body : bodies) {
        std::cout << "   "; 
        Vec<dim> pointPos = body.position;
        pointPos.print(); 
        std::cout << std::endl;
    }
    if (isLeaf) return;
    for (int childId = 0; childId < childrenCnt; childId++) {
        children[childId]->printInfoRecursively();
    }
}

template <int dim>
void BHTree<dim>::Node::calculateCentersOfMassRecursively() {
    calculateCenterOfMass();
    if (!isLeaf) {
        for (int childId = 0; childId < childrenCnt; childId++) {
            children[childId]->calculateCentersOfMassRecursively();
        }
    }
}

template <int dim>
void BHTree<dim>::Node::calculateCenterOfMass() {
    if (bodies.empty()) {
        return;
    }

    totalMass = 0.;
    centerOfMass = 0.;

    for (const Body<dim> &body : bodies) {
        totalMass += body.mass;
        centerOfMass += Vec<dim>(body.position) * body.mass;
    }

    centerOfMass /= totalMass;
}


template <int dim>
Vec<dim> BHTree<dim>::calculateForce(const Body<dim>& b) const {
    return root->calculateForce(b);
}

template <int dim>
Vec<dim> BHTree<dim>::Node::calculateForce(const Body<dim>& b) const {
    
}

// -- Explicit template instantiations --

// Instantiate for 2D
template class BHTree<2>;
template class Box<2>;

// Instantiate for 3D
template class BHTree<3>;
template class Box<3>;