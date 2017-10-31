//
//  fluidparticles.hpp
//  pbrt
//
//  Created by Thomas Roughton on 31/10/17.
//

#ifndef fluidparticles_hpp
#define fluidparticles_hpp

#include "primitive.h"

namespace pbrt {
    
    struct PositionKeyframe {
        PositionKeyframe(Point3f position, Float time) : position(position), time(time) { }
        
        Point3f position;
        Float time;
    };
    
    // MovingPrimitive Declarations
    class MovingPrimitive : public Primitive {
    public:
        // MovingPrimitive Public Methods
        MovingPrimitive(std::shared_ptr<Primitive> &primitive)
        : primitive(primitive) {}
        bool Intersect(const Ray &r, SurfaceInteraction *in) const;
        bool IntersectP(const Ray &r) const;
        bool IsProxy() const;
        Point3f Interpolate(Float t) const;
        const AreaLight *GetAreaLight() const { return nullptr; }
        const Material *GetMaterial() const { return nullptr; }
        void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                        MemoryArena &arena, TransportMode mode,
                                        bool allowMultipleLobes) const {
            LOG(FATAL) <<
            "TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
            "called";
        }
        Bounds3f WorldBound(Float startTime, Float endTime) const;
        
        std::vector<PositionKeyframe> keyframes;
        
    private:
        // MovingPrimitive Private Data
        std::shared_ptr<Primitive> primitive;
    };

    class BVHAccel;
    
class FluidContainer : public Primitive {
public:
    // TransformedPrimitive Public Methods
    FluidContainer(std::vector<std::shared_ptr<Primitive>> &primitives)
    : primitives(primitives) {}
    bool Intersect(const Ray &r, SurfaceInteraction *in) const;
    bool IntersectP(const Ray &r) const;
    bool IsProxy() const;
    const AreaLight *GetAreaLight() const { return nullptr; }
    const Material *GetMaterial() const { return nullptr; }
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const {
        LOG(FATAL) <<
        "FluidContainer::ComputeScatteringFunctions() shouldn't be "
        "called";
    }
    Bounds3f WorldBound(Float startTime, Float endTime) const;
    
private:
    // TransformedPrimitive Private Data
    std::vector<std::shared_ptr<Primitive>> primitives;
    mutable std::shared_ptr<BVHAccel> frameBVH;
    mutable Float bvhStartTime = Infinity, bvhEndTime = Infinity;
};
    
    std::shared_ptr<FluidContainer> CreateFluidContainer(const std::shared_ptr<Material> &material,
                                                         const MediumInterface &mediumInterface,
                                                         const bool isProxy,
                                                         int nParticles, int newParticlesPerFrame, float *positions, Float radius,
                                                         int nFrames, Float startFrame, Float frameStep
                                                         );
}

#endif /* fluidparticles_hpp */
