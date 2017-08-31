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
    };
    
    OctreeAccel::OctreeAccel(const std::vector<std::shared_ptr<Primitive>> &primitives) {
        // Initialize _primitiveInfo_ array for primitives
        std::vector<OctreePrimitiveInfo> primitiveInfo(primitives.size());
        for (size_t i = 0; i < primitives.size(); ++i) {
            primitiveInfo[i] = OctreePrimitiveInfo(i, primitives[i]->WorldBound());
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
    
    bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    }
    
    bool OctreeAccel::IntersectP(const Ray &ray) const {
    }
    
    std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
                                                         const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps) {
        return std::make_shared<OctreeAccel>(prims);
    }
    
}  // namespace pbrt
