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
        
        // OctreeAccel Public Methods
        OctreeAccel(const std::vector<std::shared_ptr<Primitive>> &p);
        Bounds3f WorldBound() const;
        ~OctreeAccel();
        bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
        bool IntersectP(const Ray &ray) const;
        
    private:
        // BVHAccel Private Methods
        
        int addNodePrimitives(std::shared_ptr<Primitive>* prims, const size_t primCount);
        
        Point3f computeCentroid(const OctreePrimitiveInfo* prims, const size_t primCount);
        
        bool isInChild(Point3f point, Point3f centroid, OctreeChild child);
        
        int childIndexBufferOffset(size_t base, size_t numChildren);
        
        size_t nextNodeIndex();
        
        size_t buildRecursive(OctreePrimitiveInfo* primInfos, std::shared_ptr<Primitive>* prims, const size_t primCount);
        
        // BVHAccel Private Data
        std::vector<std::shared_ptr<Primitive>> primitives;
        std::vector<OctreeNode> nodesBuffer;
        std::vector<int> childrenIndexBuffer;
    };
    
    std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
                                                   const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps);
    
}  // namespace pbrt


#endif // PBRT_ACCELERATORS_OCTREE_H
