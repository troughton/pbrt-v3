
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LIGHTS_LIGHTPROBE_H
#define PBRT_LIGHTS_LIGHTPROBE_H

// lights/lightprobe.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

namespace pbrt {

// LightProbe Declarations
class LightProbe {
  public:
    struct InfluenceArea {
        Vector3f ellipsoidInvABCSq;
        Float innerRadius;
        Float outerRadius;
        
        InfluenceArea(const Vector3f& ellipsoidABC, Float innerRadius, Float outerRadius) :
            innerRadius(innerRadius),
            outerRadius(outerRadius) {
                
                this->ellipsoidInvABCSq = Vector3f(1.0 / (ellipsoidABC.x * ellipsoidABC.x),
                                                   1.0 / (ellipsoidABC.y * ellipsoidABC.y),
                                                   1.0 / (ellipsoidABC.z * ellipsoidABC.z));
                
            }
        
        inline Float evaluateRadiusSq(const Vector3f& point) const {
            return point.x * point.x * ellipsoidInvABCSq.x +
                point.y * point.y * ellipsoidInvABCSq.y +
                point.z * point.z * ellipsoidInvABCSq.z;
        }
        
        inline Float evaluateRadius(const Vector3f& point) const {
            return sqrt(this->evaluateRadiusSq(point));
        }
    };
    
    // LightProbe Public Methods
    LightProbe(const Transform &LightToWorld, const Spectrum &L, const std::string &texmap, const Bounds3f& proxyVolume, const InfluenceArea& influenceArea);
    
    Spectrum Power(Float worldRadius) const;
    Spectrum Le(const Ray &ray) const;
    Float Pdf_Li(const Interaction &, const Vector3f &w) const;
    Float ComputeInfluenceWeight(const Point3f& point) const;

    const Bounds3f proxyVolume;
    const Point3f worldSpacePosition;
    const InfluenceArea influenceArea;
    
    const Transform LightToWorld, WorldToLight;
    std::unique_ptr<Distribution2D> distribution;
    std::unique_ptr<MIPMap<RGBSpectrum>> Lmap;
    
  private:
    // LightProbe Private Data

    Point3f worldCenter;
    Float worldRadius;
};

std::shared_ptr<LightProbe> CreateLightProbe(
    const Transform &light2world, const ParamSet &paramSet);

}  // namespace pbrt

#endif  // PBRT_LIGHTS_LIGHTPROBE_H
