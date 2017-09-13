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
        OctreeAccel(const std::vector<std::shared_ptr<Primitive>> &p, size_t depthLimit = 8, size_t maxPrimsPerNode = 8);
        Bounds3f WorldBound() const;
        ~OctreeAccel();
        bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
        bool IntersectP(const Ray &ray) const;
        
    private:
        // OctreeAccel Private Methods
        
        size_t addPrimitives(const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount);
        size_t nextNodeIndex();
        void makeLeaf(OctreeNode* node, const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount);
        size_t buildRecursive(const Bounds3f bounds, const std::vector<std::shared_ptr<Primitive>>& prims, OctreePrimitiveInfo* primInfos, const size_t primCount, size_t depth);
        
        
        bool processSubtreeP(size_t nodeIndex, const Ray& ray, const Vector3f &invDir, const int dirIsNeg[3]) const;
        
        bool processSubtree(size_t nodeIndex, const Ray& ray, const Vector3f &invDir, const int dirIsNeg[3], SurfaceInteraction *isect) const;
        
        bool traverseOctree(const Ray& ray, SurfaceInteraction *isect) const;
        inline const OctreeNode* nodeAt(size_t index) const;
        
        // OctreeAccel Private Data
        std::vector<std::shared_ptr<Primitive>> primitives;
        std::vector<OctreeNode> nodesBuffer;
        
        Bounds3f worldBound;
        Point3f centre;
        size_t depthLimit;
        size_t maxPrimsPerNode;
    };
    
    std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
                                                   const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps);
    
}  // namespace pbrt


#endif // PBRT_ACCELERATORS_OCTREE_H
