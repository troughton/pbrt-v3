//
//  octree.hpp
//  bsdftest
//
//  Created by Thomas Roughton on 31/08/17.
//

#ifndef PBRT_ACCELERATORS_OCTREE_H
#define PBRT_ACCELERATORS_OCTREE_H

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

// accelerators/bvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <atomic>

namespace pbrt {
    struct OctreeNode;
    
    // OctreeAccel Forward Declarations
    struct OctreePrimitiveInfo;
    
    // OctreeAccel Declarations
    class OctreeAccel : public Aggregate {
    public:
        // OctreeAccel Public Types
        
        enum class OctreeChild : uint8_t {
            // PosX means bit 2 is set
            // PosY means bit 1 is set
            // PosZ means bit 0 is set
            NegXNegYNegZ = 0,
            NegXNegYPosZ = 1,
            NegXPosYNegZ = 2,
            NegXPosYPosZ = 3,
            PosXNegYNegZ = 4,
            PosXNegYPosZ = 5,
            PosXPosYNegZ = 6,
            PosXPosYPosZ = 7
        };
        
        enum OctreeChildMask : uint8_t {
            // PosX means bit 2 is set
            // PosY means bit 1 is set
            // PosZ means bit 0 is set
            PosX = 0b100,
            PosY = 0b010,
            PosZ = 0b001
        };
        
        // OctreeAccel Public Methods
        OctreeAccel(const std::vector<std::shared_ptr<Primitive>> &p);
        Bounds3f WorldBound() const;
        ~OctreeAccel();
        bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
        bool IntersectP(const Ray &ray) const;
        
    private:
        // BVHAccel Private Methods
        
        size_t addPrimitives(const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount);
        Point3f computeCentroid(const OctreePrimitiveInfo* prims, const size_t primCount);
        int childIndexBufferOffset(size_t base, size_t numChildren);
        size_t nextNodeIndex();
        size_t buildRecursive(const Bounds3f bounds, const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount, size_t depth);
        
        size_t childNodeAtIndex(size_t nodeIndex, const OctreeNode &node, uint8_t childIndex) const;
        
        struct TraversalContext;
        
        size_t firstIntersectedNode(Vector3f t0, Vector3f tm) const;
        size_t newNode(Vector3f tm, size_t x, size_t y, size_t z) const;
        bool processSubtree(Vector3f t0, Vector3f t1, size_t nodeIndex, const TraversalContext& traversalContext) const;
        bool traverseOctree(const Ray& ray, SurfaceInteraction *isect) const;
        inline const OctreeNode* nodeAt(size_t index) const;
        
        // BVHAccel Private Data
        std::vector<std::shared_ptr<Primitive>> primitives;
        std::vector<OctreeNode> nodesBuffer;
        std::vector<uint32_t> childrenIndexBuffer;
        
        Bounds3f worldBound;
        Point3f centre;
        size_t depthLimit = 20;
        size_t maxPrimsPerNode = 4;
    };
    
    std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
                                                   const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps);
    
}  // namespace pbrt


#endif // PBRT_ACCELERATORS_OCTREE_H
