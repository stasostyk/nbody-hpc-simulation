#include "BHTree.hpp"

template<int dim, typename Attributes>
void BHTree<dim, Attributes>::printInfo() const {
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

template<int dim, typename Attributes>
void BHTree<dim, Attributes>::Node::insertBody(const Body<dim> &body) {
    if (isEmptyNode) {
        isEmptyNode = false;
        centerOfMass.combineWith(body);
        return;
    }

    centerOfMass.combineWith(body);
    insertToChild(body);
}

template <int dim, typename Attributes>
int BHTree<dim, Attributes>::Node::getMaxDepth() {
    int dep = depth;

    if (!isLeaf) {
        for (int childId = 0; childId < childrenCnt; childId++) {
            dep = std::max(dep, 1 + children[childId]->getMaxDepth());
        }
    }

    return dep;
}

template<int dim, typename Attributes>
void BHTree<dim, Attributes>::Node::initChild(int childId) {
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

template<int dim, typename Attributes>
void BHTree<dim, Attributes>::Node::insertToChild(const Body<dim> &body) {
    int childId = getChildIndex(body.position);

    if (!children[childId]) {
        initChild(childId);
    }

    children[childId]->insertBody(body);
}

template<int dim, typename Attributes>
void BHTree<dim, Attributes>::Node::printInfoRecursively() const {
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

template <int dim, typename Attributes>
Vec<dim> BHTree<dim, Attributes>::calculateForce(const Body<dim>& body, double theta, const forces::force<dim, Attributes> &force) const {
    return root->calculateForce(body, theta, force);
}

template <int dim, typename Attributes>
Vec<dim> BHTree<dim, Attributes>::Node::calculateForce(const Body<dim>& body, double theta, const forces::force<dim, Attributes> &force) const {

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

template <int dim, typename Attributes>
int BHTree<dim, Attributes>::Node::getChildIndex(const Vec<dim> &position) const {
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
template class BHTree<2, EmptyAttributes>;
template class Box<2>;

// Instantiate for 3D
template class BHTree<3, EmptyAttributes>;
template class Box<3>;