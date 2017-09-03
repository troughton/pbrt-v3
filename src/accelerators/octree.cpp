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
        
        this->primitives.shrink_to_fit();
        this->nodesBuffer.shrink_to_fit();
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
        
        Bounds3f constrainedBounds = pbrt::Intersect(primitiveBounds, bounds);
        
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
    
    
    bool OctreeAccel::processSubtree(size_t nodeIndex, const Ray& ray, const Vector3f &invDir,
                                     const int dirIsNeg[3], SurfaceInteraction *isect) const {
        const OctreeNode *node = this->nodeAt(nodeIndex);
        
        bool innerHit = false;
        
        // Check any primitives within the node.
        
        for (size_t i = 0; i < node->primitiveCount; i += 1) {
            size_t primitiveIndex = node->primitiveOffset + i;
            
            const std::shared_ptr<Primitive>& primitive = this->primitives[primitiveIndex];
            if (isect == nullptr) {
                innerHit |= primitive->IntersectP(ray);
            } else {
                innerHit |= primitive->Intersect(ray, isect);
            }
        }
        
        bool childHit = false;
        
        if (node->presentChildren != 0) { // Traverse its children.
            
            uint8_t remainingChildren = node->presentChildren;
            float childTMins[8];
            float childTMaxs[8];
            
            for (size_t childIndex = 0; childIndex < 8; childIndex += 1) {
                size_t childNode = node->childIndices[childIndex] + nodeIndex;
                bool childIsPresent = childNode != nodeIndex;
                if (childIsPresent) {
                    bool intersects = this->nodeAt(childNode)->bounds.IntersectP(ray, invDir, dirIsNeg, &childTMins[childIndex], &childTMaxs[childIndex]);
                    remainingChildren = (remainingChildren & ~(1 << childIndex)) | ((intersects ? 1 : 0) << childIndex); // zero out any children that don't intersect.
                }
            }
            
            float bestTMax = ray.tMax;
            while (remainingChildren != 0) {
                
                float bestTMin = std::numeric_limits<float>::infinity();
                size_t bestChildIndex = 0;
                
                for (size_t childIndex = 0; childIndex < 8; childIndex += 1) {
                    if ((remainingChildren & (1 << childIndex)) == 0) { continue; }
                    if (childTMins[childIndex] < bestTMin) {
                        bestTMin = childTMins[childIndex];
                        bestChildIndex = childIndex;
                    }
                }
                
                // Process that child.
                remainingChildren &= ~(1 << bestChildIndex);
                if (childTMins[bestChildIndex] < bestTMax) {
                    size_t childNode = node->childIndices[bestChildIndex] + nodeIndex;
                    if (this->processSubtree(childNode, ray, invDir, dirIsNeg, isect)) {
                        childHit = true;
                        bestTMax = std::min(bestTMax, ray.tMax);
                    }
                }
            }
        }
        
        return childHit || innerHit;
    }
    
    bool OctreeAccel::traverseOctree(const Ray& ray, SurfaceInteraction *isect) const {
        if (this->nodesBuffer.empty()) { return false; }
        
        float ta, tb;
        if (!this->worldBound.IntersectP(ray, &ta, &tb)) {
            return false;
        }
        
        Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
        
        return this->processSubtree(0, ray, invDir, dirIsNeg, isect);
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
