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
    // std::cout << " Inserting body " << body.bodyId << ", node dep=" << depth << std::endl;

    // if (bodies.empty()) {
    if (isEmptyNode) {
        isEmptyNode = false;
        // bodies.push_back(body);
        centerOfMass.combineWith(body);
        return;
    }

    if (isLeaf) {
        splitNode(); // Make sure there are children.
    }

    // bodies.push_back(body);
    centerOfMass.combineWith(body);
    insertToChild(body);
}

template <int dim>
int BHTree<dim>::Node::getMaxDepth() {
    int dep = depth;

    if (!isLeaf) {
        for (int childId = 0; childId < childrenCnt; childId++) {
            dep = std::max(dep, 1 + children[childId]->getMaxDepth());
        }
    }

    return dep;
}

template<int dim>
void BHTree<dim>::Node::initChild(int childId) {
    Vec<dim> corner1;
    Vec<dim> corner2;

    const Vec<dim> midPoint = (bounds.corner1 + bounds.corner2) * 0.5;

    for (int dimId = 0; dimId < dim; dimId++) {
        if ((childId & (1 << dimId)) == 0) {
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
    children[childId] = std::make_unique<Node>(childBounds, depth + 1);
    activeChildren.push_back(childId);  
}

template<int dim>
void BHTree<dim>::Node::insertToChild(const Body<dim> &body) {
    int childId = getChildIndex(body.position);

    if (!children[childId]) {
        initChild(childId);
    }

    children[childId]->insertBody(body);

    // // Find to which child node the new body belongs
    // for (int childId = 0; childId < childrenCnt; childId++) {
    //     Vec<dim> bodyPos(body.position);
    //     if (children[childId]->bounds.isPointInside(bodyPos)) {
    //         children[childId]->insertBody(body);
    //         return;
    //     }
    // }
}

template<int dim>
void BHTree<dim>::Node::splitNode() {
    if (!isLeaf) {
        // Already split.
        return;
    }
    isLeaf = false;

    // for (int bitmask = 0; bitmask < (int)childrenCnt; bitmask++) {
    //     initChild(bitmask);
    // }

    // for (const Body<dim> &body : bodies) {
        // insertToChild(body);
    // }
    insertToChild(centerOfMass); // creates a new node there
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
    std::cout << " mass: " << centerOfMass.mass << ", centerOfMass: ";
    centerOfMass.position.print();
    std::cout << std::endl;

    // std::cout << " bodies (positions): " << std::endl;
    // for (const Body<dim> &body : bodies) {
    //     std::cout << "   "; 
    //     Vec<dim> pointPos = body.position;
    //     pointPos.print(); 
    //     std::cout << std::endl;
    // }
    if (isLeaf) return;
    for (int childId = 0; childId < childrenCnt; childId++) {
        children[childId]->printInfoRecursively();
    }
}

// template <int dim>
// void BHTree<dim>::Node::calculateCentersOfMassRecursively() {
//     calculateCenterOfMass();
//     if (!isLeaf) {
//         for (int childId = 0; childId < childrenCnt; childId++) {
//             children[childId]->calculateCentersOfMassRecursively();
//         }
//     }
// }

// template <int dim>
// void BHTree<dim>::Node::calculateCenterOfMass() {
//     if (bodies.empty()) {
//         return;
//     }

//     totalMass = 0.;
//     centerOfMass = 0.;

//     for (const Body<dim> &body : bodies) {
//         totalMass += body.mass;
//         centerOfMass += Vec<dim>(body.position) * body.mass;
//     }

//     centerOfMass /= totalMass;
// }


template <int dim>
Vec<dim> BHTree<dim>::calculateForce(const Body<dim>& body, double theta, const forces::force<dim> &force) const {
    return root->calculateForce(body, theta, force);
}

template <int dim>
Vec<dim> BHTree<dim>::Node::calculateForce(const Body<dim>& body, double theta, const forces::force<dim> &force) const {

    if (isLeaf) {
        // Vec<dim> totalForce = 0.;
        
        // for (const Body<dim> &other : bodies) {
        //     if (other.bodyId == body.bodyId) {
        //         continue;
        //     }

        //     totalForce += force(body.position, body.mass, other.position, other.mass);
        // }

        // in the leaf node there should be max one body, so center of mass is the body
        if (!isEmptyNode && centerOfMass.bodyId != body.bodyId) {
            return force(body.position, body.mass, centerOfMass.position, centerOfMass.mass);
        }

        return Vec<dim>(0.);
    }

    // if (bodies.empty()) {
    if (isEmptyNode) {
        return Vec<dim>(0.);
    }

    const double distanceSqr = (centerOfMass.position - body.position).normSquared();

    const double bhVal = lengthSqr / (distanceSqr + distanceEps*distanceEps);
    // const double bhVal = lengthSqr / (distanceSqr );
    // std::cout << " body " << body.bodyId << ", tree dep=" << depth << " bhVal = " << bhVal << std::endl;
    // std::cout << "   len: " << sqrt(lengthSqr) << "  dist: " << sqrt(distanceSqr) << std::endl;

    // std::cout << "theta " << theta << std::endl;
    // Barnes-Hut criterion
    if (bhVal < theta) {
        // std::cout << " body " << body.bodyId << ", tree dep=" << depth << " bhVal = " << bhVal << std::endl;
        // std::cout << "   len: " << sqrt(lengthSqr) << "  dist: " << sqrt(distanceSqr) << std::endl;


        // Approximating far away bodies as one with computed center of mass.
        return force(body.position, body.mass, centerOfMass.position, centerOfMass.mass);
    } 

    Vec<dim> totalForce(0.);
    for (int childId : activeChildren) {
    // for (int childId = 0; childId < childrenCnt; childId++) {
        totalForce += children[childId]->calculateForce(body, theta, force);
    }

    return totalForce;
}

template <int dim> 
int BHTree<dim>::Node::getChildIndex(const Vec<dim> &position) const {
    const Vec<dim> midPoint = (bounds.corner1 + bounds.corner2) * 0.5;
    int index = 0;
    for (int dimId = 0; dimId < dim; dimId++) {
        if (position[dimId] >= midPoint[dimId] - boundCheckEps) {
            index |= (1 << dimId);
        }
    }
    return index;
}

// -- Explicit template instantiations --

// Instantiate for 2D
template class BHTree<2>;
template class Box<2>;

// Instantiate for 3D
template class BHTree<3>;
template class Box<3>;