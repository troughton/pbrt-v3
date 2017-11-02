//
//  fluidparticles.cpp
//  pbrt
//
//  Created by Thomas Roughton on 31/10/17.
//

#include <fstream>
#include <iterator>
#include "shapes/fluidparticles.h"
#include "shapes/sphere.h"
#include "accelerators/bvh.h"
#include "stats.h"
#include "efloat.h"
#include "sampling.h"
#include "paramset.h"

namespace pbrt {
    
    const Transform identityTransform;
    
    // Sphere Declarations
    class SimpleSphere : public Shape {
    public:
        // Sphere Public Methods
        SimpleSphere(Float radius)
        : Shape(&identityTransform, &identityTransform, false),
        radius(radius) {}
        Bounds3f ObjectBound() const;
        bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                       bool testAlphaTexture) const;
        bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
        Float Area() const;
        Interaction Sample(const Point2f &u, Float *pdf) const;
        Interaction Sample(const Interaction &ref, const Point2f &u,
                           Float *pdf) const;
        Float Pdf(const Interaction &ref, const Vector3f &wi) const;
        Float SolidAngle(const Point3f &p, int nSamples) const;
        
    private:
        // Sphere Private Data
        const Float radius;
    };
    
    // SimpleSphere Method Definitions
    Bounds3f SimpleSphere::ObjectBound() const {
        return Bounds3f(Point3f(-radius, -radius, -radius),
                        Point3f(radius, radius, radius));
    }
    
    bool SimpleSphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                                 bool testAlphaTexture) const {
        ProfilePhase p(Prof::ShapeIntersect);
        Float phi;
        Point3f pHit;
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        const Ray& ray = r;
        
        // Compute quadratic sphere coefficients
        
        // Initialize _EFloat_ ray coordinate values
        EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
        EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
        
        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) return false;
        
        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0) {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > ray.tMax) return false;
        }
        
        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);
        
        // Check to make sure this is a valid intersection given the number of fluid particles encountered.
        
        if (Dot(Vector3f(pHit), r.d) > 0) {
            // We're leaving the particle.
            
            r.insideFluidParticleCount -= 1;
            
        } else {
            // We're entering the particle.
            r.insideFluidParticleCount += 1;
        }
        
        
        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;
        
        const Float phiMax = 2 * Pi;
        const Float thetaMin = 0.f;
        const Float thetaMax = Pi;
        
        // Find parametric representation of sphere hit
        Float u = phi / phiMax;
        Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
        Float v = (theta - thetaMin) / (thetaMax - thetaMin);
        
        // Compute sphere $\dpdu$ and $\dpdv$
        Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
        Float invZRadius = 1 / zRadius;
        Float cosPhi = pHit.x * invZRadius;
        Float sinPhi = pHit.y * invZRadius;
        Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
        Vector3f dpdv =
        (thetaMax - thetaMin) *
        Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));
        
        // Compute sphere $\dndu$ and $\dndv$
        Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
        Vector3f d2Pduv =
        (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
        Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
        Vector3f(pHit.x, pHit.y, pHit.z);
        
        // Compute coefficients for fundamental forms
        Float E = Dot(dpdu, dpdu);
        Float F = Dot(dpdu, dpdv);
        Float G = Dot(dpdv, dpdv);
        Vector3f N = Normalize(Cross(dpdu, dpdv));
        Float e = Dot(N, d2Pduu);
        Float f = Dot(N, d2Pduv);
        Float g = Dot(N, d2Pdvv);
        
        // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
        Float invEGF2 = 1 / (E * G - F * F);
        Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
                                 (e * F - f * E) * invEGF2 * dpdv);
        Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
                                 (f * F - g * E) * invEGF2 * dpdv);
        
        // Compute error bounds for sphere intersection
        Vector3f pError = gamma(5) * Abs((Vector3f)pHit);
        
        // Initialize _SurfaceInteraction_ from parametric information
        *isect = SurfaceInteraction(pHit, pError, Point2f(u, v),
                                    -ray.d, dpdu, dpdv, dndu, dndv,
                                    ray.time, this);
        
        // Update _tHit_ for quadric intersection
        *tHit = (Float)tShapeHit;
        return true;
    }
    
    bool SimpleSphere::IntersectP(const Ray &r, bool testAlphaTexture) const {
        ProfilePhase p(Prof::ShapeIntersectP);
        Point3f pHit;
        // Transform _Ray_ to object space
        const Ray& ray = r;
        
        // Compute quadratic sphere coefficients
        
        // Initialize _EFloat_ ray coordinate values
        EFloat ox(ray.o.x, 0.f), oy(ray.o.y, 0.f), oz(ray.o.z, 0.f);
        EFloat dx(ray.d.x, 0.f), dy(ray.d.y, 0.f), dz(ray.d.z, 0.f);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
        
        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) return false;
        
        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0) {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > ray.tMax) return false;
        }
        
        return true;
    }
    
    Float SimpleSphere::Area() const { return 2 * Pi * radius; }
    
    Interaction SimpleSphere::Sample(const Point2f &u, Float *pdf) const {
        Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
        Interaction it;
        it.n = Normalize(Normal3f(pObj.x, pObj.y, pObj.z));
        if (reverseOrientation) it.n *= -1;
        // Reproject _pObj_ to sphere surface and compute _pObjError_
        pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
        it.p = pObj;
        *pdf = 1 / Area();
        return it;
    }
    
    Interaction SimpleSphere::Sample(const Interaction &ref, const Point2f &u,
                                     Float *pdf) const {
        Point3f pCenter = Point3f(0, 0, 0);
        
        // Sample uniformly on sphere if $\pt{}$ is inside it
        Point3f pOrigin =
        OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius) {
            Interaction intr = Sample(u, pdf);
            Vector3f wi = intr.p - ref.p;
            if (wi.LengthSquared() == 0)
                *pdf = 0;
            else {
                // Convert from area measure returned by Sample() call above to
                // solid angle measure.
                wi = Normalize(wi);
                *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
            }
            if (std::isinf(*pdf)) *pdf = 0.f;
            return intr;
        }
        
        // Compute coordinate system for sphere sampling
        Vector3f wc = Normalize(pCenter - ref.p);
        Vector3f wcX, wcY;
        CoordinateSystem(wc, &wcX, &wcY);
        
        // Sample sphere uniformly inside subtended cone
        
        // Compute $\theta$ and $\phi$ values for sample in cone
        Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
        Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
        Float cosTheta = (1 - u[0]) + u[0] * cosThetaMax;
        Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
        Float phi = u[1] * 2 * Pi;
        
        // Compute angle $\alpha$ from center of sphere to sampled point on surface
        Float dc = Distance(ref.p, pCenter);
        Float ds = dc * cosTheta -
        std::sqrt(std::max(
                           (Float)0, radius * radius - dc * dc * sinTheta * sinTheta));
        Float cosAlpha = (dc * dc + radius * radius - ds * ds) / (2 * dc * radius);
        Float sinAlpha = std::sqrt(std::max((Float)0, 1 - cosAlpha * cosAlpha));
        
        // Compute surface normal and sampled point on sphere
        Vector3f nWorld =
        SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
        Point3f pWorld = pCenter + radius * Point3f(nWorld.x, nWorld.y, nWorld.z);
        
        // Return _Interaction_ for sampled point on sphere
        Interaction it;
        it.p = pWorld;
        it.pError = gamma(5) * Abs((Vector3f)pWorld);
        it.n = Normal3f(nWorld);
        if (reverseOrientation) it.n *= -1;
        
        // Uniform cone PDF.
        *pdf = 1 / (2 * Pi * (1 - cosThetaMax));
        
        return it;
    }
    
    Float SimpleSphere::Pdf(const Interaction &ref, const Vector3f &wi) const {
        Point3f pCenter = Point3f(0, 0, 0);
        // Return uniform PDF if point is inside sphere
        Point3f pOrigin =
        OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
            return Shape::Pdf(ref, wi);
        
        // Compute general sphere PDF
        Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
        Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
        return UniformConePdf(cosThetaMax);
    }
    
    Float SimpleSphere::SolidAngle(const Point3f &p, int nSamples) const {
        Point3f pCenter = Point3f(0, 0, 0);
        if (DistanceSquared(p, pCenter) <= radius * radius)
            return 4 * Pi;
        Float sinTheta2 = radius * radius / DistanceSquared(p, pCenter);
        Float cosTheta = std::sqrt(std::max((Float)0, 1 - sinTheta2));
        return (2 * Pi * (1 - cosTheta));
    }
    
    
    /// MovingPrimitive
    
    Point3f MovingPrimitive::Interpolate(Float time) const {
        if (keyframes.size() == 1 || time <= keyframes.front().time) {
            return keyframes.front().position;
        }
        if (time >= keyframes.back().time) {
            return keyframes.back().position;
        }
        
        auto endKeyframeIt = std::lower_bound(keyframes.begin(), keyframes.end(), time, [](const PositionKeyframe &lhs, Float rhs) {
            return lhs.time < rhs;
        });
        auto beginKeyframeIt = (endKeyframeIt == keyframes.cbegin()) ? endKeyframeIt : (endKeyframeIt - 1);
        
        Float dt = (time - beginKeyframeIt->time) / (endKeyframeIt->time - beginKeyframeIt->time);
        
        // Interpolate translation at _dt_
        return (1 - dt) * beginKeyframeIt->position + dt * endKeyframeIt->position;
    }
    
    Bounds3f MovingPrimitive::WorldBound(Float startTime, Float endTime) const {
        Bounds3f primitiveBound = primitive->WorldBound(startTime, endTime);
        
        auto beginKeyframeIt = std::lower_bound(keyframes.begin(), keyframes.end(), startTime, [](const PositionKeyframe &lhs, Float rhs) {
            return lhs.time < rhs;
        });
        
        if (beginKeyframeIt != keyframes.begin()) {
            beginKeyframeIt--;
        }
        
        auto endKeyframeIt = std::lower_bound(keyframes.begin(), keyframes.end(), endTime, [](const PositionKeyframe &lhs, Float rhs) {
            return lhs.time < rhs;
        });
        
        if (endKeyframeIt == keyframes.end()) {
            endKeyframeIt--;
        }
        
        Bounds3f bounds(primitiveBound.pMin + beginKeyframeIt->position, primitiveBound.pMax + beginKeyframeIt->position);
        
        while (beginKeyframeIt < endKeyframeIt) {
            auto secondKeyframeIt = beginKeyframeIt + 1;
            bounds = Union(bounds, Bounds3f(primitiveBound.pMin + secondKeyframeIt->position, primitiveBound.pMax + secondKeyframeIt->position));
            
            beginKeyframeIt++;
        }
        return bounds;
    }
    
    bool MovingPrimitive::IsProxy() const {
        return primitive->IsProxy();
    }
    
    bool MovingPrimitive::IntersectP(const Ray &r) const {
        Point3f positionOffset = this->Interpolate(r.time);
        Ray ray = r;
        ray.o = ray.o + (-positionOffset);
        return primitive->IntersectP(ray);
    }
    
    bool MovingPrimitive::Intersect(const Ray &r,
                                    SurfaceInteraction *isect) const {
        
        Point3f positionOffset = this->Interpolate(r.time);
        Ray ray = r;
        ray.o = ray.o + (-positionOffset);
        
        if (!primitive->Intersect(ray, isect)) return false;
        
        r.tMax = ray.tMax;
        r.insideFluidParticleCount = ray.insideFluidParticleCount;
        
        isect->p += positionOffset;
        
        return true;
    }
    
    
    /// FluidContainer
    
    Bounds3f FluidContainer::WorldBound(Float startTime, Float endTime) const {
        if (startTime != this->bvhStartTime || endTime != this->bvhEndTime) {
            // Construct the frame BVH
            
            auto primitivesEnd = std::lower_bound(this->primitives.begin(), this->primitives.end(), startTime, [](std::shared_ptr<Primitive> prim, Float endTime) -> bool {
                return ((MovingPrimitive*)prim.get())->keyframes.front().time < endTime;
            });
            
            std::vector<std::shared_ptr<Primitive>> validPrimitives(primitives.begin(), primitivesEnd);
            
            ParamSet paramSet;
            this->frameBVH = CreateBVHAccelerator(validPrimitives, startTime, endTime, paramSet);
            this->bvhStartTime = startTime;
            this->bvhEndTime = endTime;
        }
        return this->frameBVH->WorldBound(startTime, endTime);
    }
    
    bool FluidContainer::IntersectP(const Ray &r) const {
        return this->frameBVH->IntersectP(r);
    }
    
    bool FluidContainer::Intersect(const Ray &r,
                                   SurfaceInteraction *isect) const {
        Ray ray = r;
        int inputInsideFluidParticleCount = ray.insideFluidParticleCount;
        
        if (this->frameBVH->Intersect(ray, isect)) {
            if (inputInsideFluidParticleCount != 0) {
                while (r.insideFluidParticleCount != 0) {
                    ray.o = ray(ray.tMax + MachineEpsilon);
                    ray.tMax = Infinity;
                    if (!this->frameBVH->Intersect(ray, isect)) {
                        break;
                    }
                }
            }
            
            r.tMax = ray.tMax;
            r.insideFluidParticleCount = ray.insideFluidParticleCount;
            
            return true;
        }
        
        return false;
    }
    
    bool FluidContainer::IsProxy() const {
        return this->frameBVH->IsProxy();
    }
    
    
    std::shared_ptr<FluidContainer> CreateFluidContainer(const std::shared_ptr<Material> &material,
                                                         const MediumInterface &mediumInterface,
                                                         const bool isProxy,
                                                         const ParamSet &params
                                                         ) {
        
        int nParticles = params.FindOneInt("particlecount", 0);
        int newParticlesPerFrame = params.FindOneInt("newparticlesperframe", 0);
        Float radius = params.FindOneFloat("particleradius", 0.05);
        std::string positionsFile = params.FindOneFilename("positionsfile", "");
        int nFrames = params.FindOneInt("framecount", 0);
        Float startFrame = params.FindOneFloat("starttime", 0);
        Float frameStep = params.FindOneFloat("framestep", 1.0);
        
        if (positionsFile.empty()) {
            Error("No data file specified for fluid container.");
            exit(-1);
        }
        
        
        
        
        std::shared_ptr<Shape> sphere = std::make_shared<SimpleSphere>(radius);
        std::shared_ptr<Primitive> spherePrim = std::make_shared<GeometricPrimitive>(sphere, material, nullptr, mediumInterface, isProxy);
        
        std::ifstream input( positionsFile, std::ios::binary );
        
        std::vector<char> positionsBuffer((
                                  std::istreambuf_iterator<char>(input)),
                                 (std::istreambuf_iterator<char>()));
        float *positions = (float*)positionsBuffer.data();
        
        std::vector<std::shared_ptr<Primitive>> particlePrims;
        
        for (int i = 0; i < nParticles; i += 1) {
            particlePrims.push_back(std::make_shared<MovingPrimitive>(spherePrim));
        }
        
        int numParticles = 0;
        
        for (int frame = 0; frame < nFrames; frame += 1) {
            numParticles += newParticlesPerFrame;
            numParticles = std::min(numParticles, nParticles);
            
            Float frameTime = startFrame + frame * frameStep;
            
            for (int i = 0; i < numParticles; i += 1) {
                MovingPrimitive* prim = (MovingPrimitive*)particlePrims[i].get();
                
                Point3f position((Float)positions[0], (Float)positions[1], (Float)positions[2]);
                prim->keyframes.push_back(PositionKeyframe(position, frameTime));
                positions += 3;
            }
        }
        
        return std::make_shared<FluidContainer>(particlePrims);
    }
}
