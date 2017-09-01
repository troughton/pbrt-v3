//
//  octree.cpp
//  bsdftest
//
//  Created by Thomas Roughton on 31/08/17.
//

// accelerators/octree.cpp*
#include "accelerators/octree.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>

namespace pbrt {
    
    //    STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
    //    STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
    //    STAT_COUNTER("BVH/Interior nodes", interiorNodes);
    //    STAT_COUNTER("BVH/Leaf nodes", leafNodes);
    
    // OctreeAccel Local Declarations
    struct OctreePrimitiveInfo {
        OctreePrimitiveInfo() {}
        OctreePrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
        : primitiveNumber(primitiveNumber),
        bounds(bounds),
        centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}
        size_t primitiveNumber;
        Bounds3f bounds;
        Point3f centroid;
    };
    
    struct OctreeNode {
        Point3f centroid;
        uint8_t presentChildren = 0;
        uint8_t primitiveCount = 0;
        int primitiveOffset = 0;
        int childrenIndexBufferOffset = 0;
        
        inline bool isTerminal() const {
            return this->presentChildren == 0;
        }
    };
    
    OctreeAccel::OctreeAccel(const std::vector<std::shared_ptr<Primitive>> &primitives) {
        // Initialize _primitiveInfo_ array for primitives
        
        std::vector<OctreePrimitiveInfo> primitiveInfo(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i) {
            primitiveInfo[i] = OctreePrimitiveInfo(i, primitives[i]->WorldBound());
            this->worldBound = Union(this->worldBound, primitiveInfo[i].bounds);
        }
        
        std::vector<std::shared_ptr<Primitive>> primitivesCopy = primitives;
        
        this->buildRecursive(primitiveInfo.data(), primitivesCopy.data(), primitives.size());
    }
    
    OctreeAccel::~OctreeAccel() { }
    
    int OctreeAccel::addNodePrimitives(std::shared_ptr<Primitive>* prims, const size_t primCount) {
        if (primCount == 0) { return 0; }
        
        size_t index = this->primitives.size();
        this->primitives.resize(index + primCount);
        
        memcpy(&this->primitives[index], prims, primCount * sizeof(prims[0]));
        return static_cast<int>(index);
    }
    
    Point3f OctreeAccel::computeCentroid(const OctreePrimitiveInfo* prims, const size_t primCount) {
        Point3f centroid(0, 0, 0);
        
        for (size_t i = 0; i < primCount; i += 1) {
            centroid += prims[i].centroid;
        }
        
        return centroid / primCount;
    }
    
    bool OctreeAccel::isInChild(Point3f point, Point3f centroid, OctreeChild child) {
        uint8_t mask = 0;
        mask |= (point.x >= centroid.x ? 0b100 : 0);
        mask |= (point.y >= centroid.y ? 0b010 : 0);
        mask |= (point.z >= centroid.z ? 0b001 : 0);
        
        return (OctreeChild)mask == child;
    }
    
    int OctreeAccel::childIndexBufferOffset(size_t base, size_t numChildren) {
        size_t returnValue = this->childrenIndexBuffer.size() - base;
        size_t requiredSize = this->childrenIndexBuffer.size() + numChildren;
        if (requiredSize > this->childrenIndexBuffer.capacity()) {
            this->childrenIndexBuffer.reserve(requiredSize);
        }
        this->childrenIndexBuffer.resize(requiredSize);
        
        return static_cast<int>(returnValue);
    }
    
    size_t OctreeAccel::nextNodeIndex() {
        size_t index = this->nodesBuffer.size();
        this->nodesBuffer.push_back(OctreeNode());
        return index;
    }
    
    size_t OctreeAccel::buildRecursive(OctreePrimitiveInfo* primInfos, std::shared_ptr<Primitive>* prims, const size_t primCount) {
        const size_t nodeIndex = this->nextNodeIndex();
        
        Point3f centroid = this->computeCentroid(primInfos, primCount);
        
        size_t categorisedPrimitivesDivider = 0;
        size_t rangeStartIndices[9];
        size_t childCount = 0;
        
        rangeStartIndices[8] = primCount;
        
        int *childIndexBuffer = nullptr;
        {
            OctreeNode *node = &this->nodesBuffer[nodeIndex];
            node->centroid = centroid;
            
            for (size_t i = 0; i < primCount; i += 1) {
                if (Inside(centroid, primInfos[i].bounds)) {
                    std::swap(prims[i], prims[categorisedPrimitivesDivider]); // Everything before the divider has already been categorised.
                    std::swap(primInfos[i], primInfos[categorisedPrimitivesDivider]);
                    categorisedPrimitivesDivider += 1;
                }
            }
            
            node->primitiveCount = categorisedPrimitivesDivider;
            
            // Everything in the range 0..<categorisedPrimitivesDivider belongs to node.
            node->primitiveOffset = this->addNodePrimitives(prims, categorisedPrimitivesDivider);
            
            // Now, assign each remaining primitive to one of the children.
            for (size_t childIndex = 0; childIndex < 8; childIndex += 1) {
                OctreeChild child = (OctreeChild)childIndex;
                size_t rangeStart = categorisedPrimitivesDivider;
                rangeStartIndices[childIndex] = rangeStart;
                for (size_t i = rangeStart; i < primCount; i += 1) {
                    if (this->isInChild(primInfos[i].centroid, centroid, child)) {
                        std::swap(prims[i], prims[categorisedPrimitivesDivider]); // Everything before the divider has already been categorised.
                        std::swap(primInfos[i], primInfos[categorisedPrimitivesDivider]);
                        categorisedPrimitivesDivider += 1;
                    }
                }
                
                // Everything in the range rangeStart..<categorisedPrimitivesDivider belongs to child.
                
                if (rangeStart != categorisedPrimitivesDivider) {
                    node->presentChildren |= 1 << childIndex;
                    childCount += 1;
                }
            }
            
            node->childrenIndexBufferOffset = this->childIndexBufferOffset(nodeIndex, childCount);
            childIndexBuffer = &this->childrenIndexBuffer[nodeIndex + node->childrenIndexBufferOffset];
        }
        
        // Now actually create the children and add them to the children index buffer.
        // After this point, the node pointer may be invalidated, so don't access it.
        
        size_t childNumber = 0;
        
        for (size_t childIndex = 0; childIndex < 8; childIndex += 1) {
            size_t rangeStart = rangeStartIndices[childIndex]; // inclusive
            size_t rangeUpperBound = rangeStartIndices[childIndex + 1]; // exclusive
            size_t rangeSize = rangeUpperBound - rangeStart;
            
            if (rangeSize == 0) { continue; };
            
            childIndexBuffer[childNumber] = this->buildRecursive(&primInfos[rangeStart], &prims[rangeStart], rangeSize);
            
            childNumber += 1;
        }
        
        return nodeIndex;
    }
    
    // Used to lookup the number of other child nodes before the current child node. More directly, hammingWeightTable[i] is the number of set bits in i given i < 128.
    uint8_t HammingWeightTable[128] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7 };
    
    size_t OctreeAccel::childNodeAtIndex(size_t nodeIndex, const OctreeNode &node, uint8_t childIndex) const {
        // The offset in the child indices buffer is the number of children present before childIndex.
        uint8_t mask = ~(~((uint8_t)0) << childIndex); // such that bits 0..<childIndex are 1 and all other bits are 0.
        uint8_t maskedValue = node.presentChildren & mask;
        
        // Now, we need to know how many bits are set in maskedValue.
        uint8_t priorChildren = HammingWeightTable[maskedValue];
        size_t childIndexInBuffer = nodeIndex + node.childrenIndexBufferOffset + priorChildren;
        return this->childrenIndexBuffer[childIndexInBuffer];
    }
    
    // Octree traversal adapted from https://stackoverflow.com/questions/10228690/ray-octree-intersection-algorithms
    
    struct OctreeAccel::TraversalContext {
        const Ray& originalRay;
        uint8_t childMask;
        Point3f rayOrigin;
        Vector3f rayInvDir;
        SurfaceInteraction *isect = nullptr;
        
        TraversalContext(const Ray& originalRay) : originalRay(originalRay) { }
    };
    
    size_t OctreeAccel::firstIntersectedNode(Vector3f t0, Vector3f tm) const {
        uint8_t answer = 0;   // initialize to 00000000
        // select the entry plane and set bits
        if (t0.x > t0.y) {
            if (t0.x > t0.z) { // PLANE YZ
                if (tm.y < t0.x) { answer |= 2; }    // set bit at position 1
                if (tm.z < t0.x) { answer |= 1; }    // set bit at position 0
                return static_cast<size_t>(answer);
            }
        } else {
            if (t0.y > t0.z){ // PLANE XZ
                if (tm.x < t0.y) { answer |= 4; }    // set bit at position 2
                if (tm.z < t0.y) { answer |= 1; }    // set bit at position 0
                return static_cast<size_t>(answer);
            }
        }
        // PLANE XY
        if (tm.x < t0.z) { answer |= 4; }    // set bit at position 2
        if (tm.y < t0.z) { answer |= 2; }    // set bit at position 1
        return static_cast<size_t>(answer);
    }
    
    size_t OctreeAccel::newNode(Vector3f tm, size_t x, size_t y, size_t z) const {
        if (tm.x < tm.y) {
            if (tm.x < tm.z) { return x; }  // YZ plane
        } else {
            if (tm.y < tm.z) { return y; } // XZ plane
        }
        return z; // XY plane;
    }
    
    bool OctreeAccel::processSubtree(Vector3f t0, Vector3f t1, size_t nodeIndex, const TraversalContext& traversalContext) const {
        if (t1.x < 0 || t1.y < 0 || t1.z < 0) { return false; }
        
        const OctreeNode *node = this->nodeAt(nodeIndex);
        
        Vector3f offset(node->centroid - traversalContext.rayOrigin);
        Vector3f tm(traversalContext.rayInvDir.x * offset.x, traversalContext.rayInvDir.y * offset.y, traversalContext.rayInvDir.z * offset.z);
        
        const uint8_t a = traversalContext.childMask;
        
        bool hit = false;
        
        // Process any primitives within this node.
        
        for (size_t i = 0; i < node->primitiveCount; i += 1) {
            size_t primitiveIndex = node->primitiveOffset + i;
            const std::shared_ptr<Primitive>& primitive = this->primitives[primitiveIndex];
            if (traversalContext.isect == nullptr) {
                hit |= primitive->IntersectP(traversalContext.originalRay);
            } else {
                hit |= primitive->Intersect(traversalContext.originalRay, traversalContext.isect);
            }
        }
        
        if (node->presentChildren == 0) {
            return hit;
        }
        
        size_t currentNode = this->firstIntersectedNode(t0, tm);
        
        do {
            switch (currentNode) {
                case 0: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, a);
                    hit |= this->processSubtree(t0, tm, childNode, traversalContext);
                    currentNode = this->newNode(tm, 4, 2, 1);
                    break;
                }
                case 1: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, 1 ^ a);
                    Vector3f newT1(tm.x, tm.y, t1.z);
                    hit |= this->processSubtree(Vector3f(t0.x, t0.y, tm.z), newT1, childNode, traversalContext);
                    currentNode = this->newNode(newT1, 5, 3, 8);
                    break;
                }
                case 2: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, 2 ^ a);
                    Vector3f newT1(tm.x, t1.y, tm.z);
                    hit |= this->processSubtree(Vector3f(t0.x, tm.y, t0.z), newT1, childNode, traversalContext);
                    currentNode = this->newNode(newT1, 6, 8, 3);
                    break;
                }
                case 3: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, 3 ^ a);
                    Vector3f newT1(tm.x, t1.y, t1.z);
                    hit |= this->processSubtree(Vector3f(t0.x, tm.y, tm.z), newT1, childNode, traversalContext);
                    currentNode = this->newNode(newT1, 7, 8, 8);
                    break;
                }
                case 4: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, 4 ^ a);
                    Vector3f newT1(t1.x, tm.y, tm.z);
                    hit |= this->processSubtree(Vector3f(tm.x, t0.y, t0.z), newT1, childNode, traversalContext);
                    currentNode = this->newNode(newT1, 8, 6, 5);
                    break;
                }
                case 5: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, 5 ^ a);
                    Vector3f newT1(t1.x, tm.y, t1.z);
                    hit |= this->processSubtree(Vector3f(tm.x, t0.y, tm.z), newT1, childNode, traversalContext);
                    currentNode = this->newNode(newT1, 8, 7, 8);
                    break;
                }
                case 6: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, 6 ^ a);
                    Vector3f newT1(t1.x, t1.y, tm.z);
                    hit |= this->processSubtree(Vector3f(tm.x, tm.y, t0.z), newT1, childNode, traversalContext);
                    currentNode = this->newNode(newT1, 8, 8, 7);
                    break;
                }
                case 7: {
                    size_t childNode = this->childNodeAtIndex(nodeIndex, *node, 7 ^ a);
                    hit |= this->processSubtree(tm, t1, childNode, traversalContext);
                    currentNode = 8;
                    break;
                }
            }
        } while (currentNode < 8);
        
        return hit;
    }
    
    bool OctreeAccel::traverseOctree(const Ray& ray, SurfaceInteraction *isect) const {
        if (this->nodesBuffer.empty()) { return false; }
        
        uint8_t childMask = 0;
        
        Ray tmpRay = ray;
        // fixes for rays with negative direction
        if (tmpRay.d.x < 0){
            tmpRay.o.x = this->centre.x * 2 - ray.o.x;
            tmpRay.d.x = -ray.d.x;
            childMask |= 4; // set bit x
        }
        
        if (tmpRay.d.y < 0) {
            tmpRay.o.y = this->centre.y * 2 - ray.o.y;
            tmpRay.d.y = -ray.d.y;
            childMask |= 2; // set bit y
        }
        
        if (tmpRay.d.z < 0) {
            tmpRay.o.z = this->centre.z * 2 - ray.o.z;
            tmpRay.d.z = -ray.d.z;
            childMask |= 1; // set bit z
        }
        
        Vector3f invD = Vector3f(1.0 / tmpRay.d.x, 1.0 / tmpRay.d.y, 1.0 / tmpRay.d.z);
        
        Vector3f t0 = (this->worldBound.pMin - tmpRay.o);
        Vector3f t1 = (this->worldBound.pMax - tmpRay.o);
        t0.x /= invD.x;
        t1.x /= invD.x;
        t0.y /= invD.y;
        t1.y /= invD.y;
        t0.z /= invD.z;
        t1.z /= invD.z;
        
        if (std::max(std::max(t0.x, t0.y), t0.z) < std::min(std::min(t1.x, t1.y), t1.z)) {
            TraversalContext context(ray);
            context.childMask = childMask;
            context.rayOrigin = tmpRay.o;
            context.rayInvDir = invD;
            context.isect = isect;
            
            return this->processSubtree(t0, t1, 0, context);
        }
        
        return false;
    }
    
    bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
        ProfilePhase p(Prof::AccelIntersect);
        return this->traverseOctree(ray, isect);
    }
    
    bool OctreeAccel::IntersectP(const Ray &ray) const {
        ProfilePhase p(Prof::AccelIntersectP);
        return this->traverseOctree(ray, nullptr);
    }
    
    std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
                                                         const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps) {
        return std::make_shared<OctreeAccel>(prims);
    }
    
}  // namespace pbrt
