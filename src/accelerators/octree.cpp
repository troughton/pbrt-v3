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
    
    STAT_COUNTER("Octree/Depth", maxDepth);
    STAT_RATIO("Octree/Primitives per leaf node", totalPrimitives, totalLeafNodes);
    
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
        Bounds3f bounds;
        uint32_t primitiveOffset = 0;
        uint16_t primitiveCount = 0;
        uint16_t presentChildren = 0;
        uint16_t childIndices[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        
        OctreeNode() { }
    };
    
    OctreeAccel::OctreeAccel(const std::vector<std::shared_ptr<Primitive>> &primitives) {
        ProfilePhase _(Prof::AccelConstruction);
        
        // Initialize _primitiveInfo_ array for primitives
        std::vector<OctreePrimitiveInfo> primitiveInfo(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i) {
            primitiveInfo[i] = OctreePrimitiveInfo(i, primitives[i]->WorldBound());
            this->worldBound = Union(this->worldBound, primitiveInfo[i].bounds);
        }
        
        totalPrimitives = primitives.size();
        
        this->primitives.reserve(primitives.size());
        
        this->nodesBuffer.reserve(2 * primitives.size());
        
        this->buildRecursive(this->worldBound, primitives, primitiveInfo.data(), primitiveInfo.size(), 0);
    }
    
    OctreeAccel::~OctreeAccel() { }
    
    Bounds3f OctreeAccel::WorldBound() const {
        return this->worldBound;
    }
    
    size_t OctreeAccel::addPrimitives(const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount) {
        size_t baseIndex = this->primitives.size();
        
        if (baseIndex + primCount > this->primitives.capacity()) {
            this->primitives.reserve(2 * (baseIndex + primCount));
        }
        
        for (size_t i = 0; i < primCount; i += 1) {
            size_t primitiveIndex = primInfos[i].primitiveNumber;
            this->primitives.push_back(prims[primitiveIndex]);
        }
        
        return baseIndex;
    }
    
    size_t OctreeAccel::nextNodeIndex() {
        size_t index = this->nodesBuffer.size();
        
        if (index > this->nodesBuffer.capacity()) {
            this->nodesBuffer.reserve(2 * index);
        }
        
        this->nodesBuffer.push_back(OctreeNode());
        
        return index;
    }
    
    void OctreeAccel::makeLeaf(OctreeNode* node, const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount) {
        node->presentChildren = 0;
        node->primitiveCount = primCount;
        node->primitiveOffset = this->addPrimitives(prims, primInfos, primCount);
        totalLeafNodes += 1;
    }
    
    size_t OctreeAccel::buildRecursive(const Bounds3f bounds, const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount, size_t depth) {
        maxDepth = std::max(maxDepth, (int64_t)depth);
        
        const size_t nodeIndex = this->nextNodeIndex();
        OctreeNode *node = &this->nodesBuffer[nodeIndex];
        
        Bounds3f primitiveBounds = primInfos[0].bounds;
        
        for (size_t i = 1; i < primCount; i += 1) {
            primitiveBounds = Union(primInfos[i].bounds, primitiveBounds);
        }
        
        Bounds3f constrainedBounds = bounds; // pbrt::Intersect(primitiveBounds, bounds);
        
        node->bounds = constrainedBounds;
        const Point3f centroid = (constrainedBounds.pMin + constrainedBounds.pMax) * 0.5;
        
        // If we're at the stopping criteria (e.g. primCount <= maxPrimsPerLeaf), stop and fill in the leaf node.
        // Otherwise, split this node into its children.
        // For each child, reorder the primInfos such that the start of the buffer contains only those elements that need to be considered for the current child.
        
        if (constrainedBounds.Volume() < 1e-2 || depth == depthLimit || primCount <= maxPrimsPerNode) {
            // Make this a leaf node.
            this->makeLeaf(node, prims, primInfos, primCount);
            
            return nodeIndex;
        }
        
        size_t nodePrimitivesEndIndex = 0;
        
        for (size_t i = 0; i < primCount; i += 1) {
            if (Inside(centroid, primInfos[i].bounds)) {
                std::swap(primInfos[i], primInfos[nodePrimitivesEndIndex]);
                nodePrimitivesEndIndex += 1;
            }
        }
        
        if (nodePrimitivesEndIndex > 0) {
            this->makeLeaf(node, prims, primInfos, nodePrimitivesEndIndex);
        }
        
        for (size_t child = 0; child < 8; child += 1) {
            bool posX = (child & OctreeChildMask::PosX) != 0;
            bool posY = (child & OctreeChildMask::PosY) != 0;
            bool posZ = (child & OctreeChildMask::PosZ) != 0;
            const Point3f minPoint = Point3f(posX ? centroid.x : node->bounds.pMin.x, posY ? centroid.y : node->bounds.pMin.y, posZ ? centroid.z : node->bounds.pMin.z);
            const Point3f maxPoint = Point3f(posX ? node->bounds.pMax.x : centroid.x, posY ? node->bounds.pMax.y : centroid.y, posZ ? node->bounds.pMax.z : centroid.z);
            
            const Bounds3f childBounds = Bounds3f(minPoint, maxPoint);
            
            size_t rangeEndIndex = nodePrimitivesEndIndex;
            
            for (size_t i = nodePrimitivesEndIndex; i < primCount; i += 1) {
                if (Overlaps(primInfos[i].bounds, childBounds)) {
                    std::swap(primInfos[i], primInfos[rangeEndIndex]);
                    rangeEndIndex += 1;
                }
            }
            
            size_t rangeSize = rangeEndIndex - nodePrimitivesEndIndex;
            
            if (rangeSize > 0) {
                size_t childNode = this->buildRecursive(childBounds, prims, &primInfos[nodePrimitivesEndIndex], rangeSize, depth + 1);
                node = &this->nodesBuffer[nodeIndex]; // Revalidate our node pointer in case of buffer resizing.
                node->childIndices[child] = childNode - nodeIndex;
                node->presentChildren |= 1 << child;
            }
        }
        
        return nodeIndex;
    }
    
    const OctreeNode* OctreeAccel::nodeAt(size_t index) const {
        return &this->nodesBuffer[index];
    }
    
    // Octree traversal adapted from https://stackoverflow.com/questions/10228690/ray-octree-intersection-algorithms
    
    struct OctreeAccel::TraversalContext {
        const Ray& originalRay;
        Vector3f origInvDir;
        int origDirIsNeg[3];
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
                if (tm.y < t0.x) { answer |= OctreeChildMask::PosY; }    // set bit at position 1
                if (tm.z < t0.x) { answer |= OctreeChildMask::PosZ; }    // set bit at position 0
                return static_cast<size_t>(answer);
            }
        } else {
            if (t0.y > t0.z){ // PLANE XZ
                if (tm.x < t0.y) { answer |= OctreeChildMask::PosX; }    // set bit at position 2
                if (tm.z < t0.y) { answer |= OctreeChildMask::PosZ; }    // set bit at position 0
                return static_cast<size_t>(answer);
            }
        }
        // PLANE XY
        if (tm.x < t0.z) { answer |= OctreeChildMask::PosX; }    // set bit at position 2
        if (tm.y < t0.z) { answer |= OctreeChildMask::PosY; }    // set bit at position 1
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
    
    
    bool OctreeAccel::processSubtree(Vector3f t0, Vector3f t1, size_t nodeIndex, TraversalContext& traversalContext) const {
        if (t1.x < 0 || t1.y < 0 || t1.z < 0) { return false; }
        
        const OctreeNode *node = this->nodeAt(nodeIndex);
        
        if (!node->bounds.IntersectP(traversalContext.originalRay, traversalContext.origInvDir, traversalContext.origDirIsNeg)) {
            return false;
        }
        
        bool hit = false;
        
        if (node->presentChildren != 0) { // Traverse its children.
            
            Vector3f tm = (t0 + t1) * 0.5;
            const uint8_t negationMask = traversalContext.childMask;
            
            size_t currentNode = this->firstIntersectedNode(t0, tm);
            
            do {
                uint8_t childIndex = currentNode ^ negationMask;
                size_t childNode = node->childIndices[childIndex] + nodeIndex;
                bool childIsPresent = childNode != nodeIndex; // The octree is sparse, so we don't always have every child.
                
                switch (currentNode) {
                    case 0: {
                        if (childIsPresent) { hit |= this->processSubtree(t0, tm, childNode, traversalContext); }
                        currentNode = this->newNode(tm, 4, 2, 1);
                        break;
                    }
                    case 1: {
                        Vector3f newT1(tm.x, tm.y, t1.z);
                        if (childIsPresent) { hit |= this->processSubtree(Vector3f(t0.x, t0.y, tm.z), newT1, childNode, traversalContext); }
                        currentNode = this->newNode(newT1, 5, 3, 8);
                        break;
                    }
                    case 2: {
                        Vector3f newT1(tm.x, t1.y, tm.z);
                        if (childIsPresent) { hit |= this->processSubtree(Vector3f(t0.x, tm.y, t0.z), newT1, childNode, traversalContext); }
                        currentNode = this->newNode(newT1, 6, 8, 3);
                        break;
                    }
                    case 3: {
                        Vector3f newT1(tm.x, t1.y, t1.z);
                        if (childIsPresent) { hit |= this->processSubtree(Vector3f(t0.x, tm.y, tm.z), newT1, childNode, traversalContext); }
                        currentNode = this->newNode(newT1, 7, 8, 8);
                        break;
                    }
                    case 4: {
                        Vector3f newT1(t1.x, tm.y, tm.z);
                        if (childIsPresent) { hit |= this->processSubtree(Vector3f(tm.x, t0.y, t0.z), newT1, childNode, traversalContext); }
                        currentNode = this->newNode(newT1, 8, 6, 5);
                        break;
                    }
                    case 5: {
                        Vector3f newT1(t1.x, tm.y, t1.z);
                        if (childIsPresent) { hit |= this->processSubtree(Vector3f(tm.x, t0.y, tm.z), newT1, childNode, traversalContext); }
                        currentNode = this->newNode(newT1, 8, 7, 8);
                        break;
                    }
                    case 6: {
                        Vector3f newT1(t1.x, t1.y, tm.z);
                        if (childIsPresent) { hit |= this->processSubtree(Vector3f(tm.x, tm.y, t0.z), newT1, childNode, traversalContext); }
                        currentNode = this->newNode(newT1, 8, 8, 7);
                        break;
                    }
                    case 7: {
                        if (childIsPresent) { hit |= this->processSubtree(tm, t1, childNode, traversalContext); }
                        currentNode = 8;
                        break;
                    }
                }
                
            } while (currentNode < 8);
        }
        
        // Check any primitives within the node.
        
        for (size_t i = 0; i < node->primitiveCount; i += 1) {
            size_t primitiveIndex = node->primitiveOffset + i;
            
            const std::shared_ptr<Primitive>& primitive = this->primitives[primitiveIndex];
            if (traversalContext.isect == nullptr) {
                hit |= primitive->IntersectP(traversalContext.originalRay);
            } else {
                hit |= primitive->Intersect(traversalContext.originalRay, traversalContext.isect);
            }
        }
        
        return hit;
    }
    
    bool OctreeAccel::traverseOctree(const Ray& ray, SurfaceInteraction *isect) const {
        if (this->nodesBuffer.empty()) { return false; }
        
        float ta, tb;
        if (!this->worldBound.IntersectP(ray, &ta, &tb)) {
            return false;
        }
        
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
        
        Vector3f invD = Vector3f(1.0 / std::max(tmpRay.d.x, MachineEpsilon), 1.0 / std::max(tmpRay.d.y, MachineEpsilon), 1.0 / std::max(tmpRay.d.z, MachineEpsilon)); // Add a small epsilon to prevent infinite (and later NaN) values.
        
        Vector3f t0 = (this->worldBound.pMin - tmpRay.o);
        Vector3f t1 = (this->worldBound.pMax - tmpRay.o);
        t0.x *= invD.x;
        t1.x *= invD.x;
        t0.y *= invD.y;
        t1.y *= invD.y;
        t0.z *= invD.z;
        t1.z *= invD.z;
        
        if (std::max(std::max(t0.x, t0.y), t0.z) < std::min(std::min(t1.x, t1.y), t1.z)) {
            TraversalContext context(ray);
            context.childMask = childMask;
            context.rayOrigin = tmpRay.o;
            context.rayInvDir = invD;
            context.isect = isect;
            
            context.origInvDir = Vector3f(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
            context.origDirIsNeg[0] = context.origInvDir.x < 0;
            context.origDirIsNeg[1] = context.origInvDir.y < 0;
            context.origDirIsNeg[2] = context.origInvDir.z < 0;
            
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
