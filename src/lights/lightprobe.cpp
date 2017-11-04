
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

// lights/infinite.cpp*
#include "lights/lightprobe.h"
#include "imageio.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {
    
    // LightProbe Method Definitions
    LightProbe::LightProbe(const Transform &LightToWorld, const Spectrum &L, const std::string &texmap, const Bounds3f& proxyVolume, const InfluenceArea &influenceArea) : proxyVolume(proxyVolume), worldSpacePosition(LightToWorld(Point3f())), influenceArea(influenceArea) {
        // Read texel data from _texmap_ and initialize _Lmap_
        Point2i resolution;
        std::unique_ptr<RGBSpectrum[]> texels(nullptr);
        if (texmap != "") {
            texels = ReadImage(texmap, &resolution);
            if (texels)
                for (int i = 0; i < resolution.x * resolution.y; ++i)
                    texels[i] *= L.ToRGBSpectrum();
        }
        if (!texels) {
            resolution.x = resolution.y = 1;
            texels = std::unique_ptr<RGBSpectrum[]>(new RGBSpectrum[1]);
            texels[0] = L.ToRGBSpectrum();
        }
        Lmap.reset(new MIPMap<RGBSpectrum>(resolution, texels.get()));
        
        // Initialize sampling PDFs for infinite area light
        
        // Compute scalar-valued image _img_ from environment map
        int width = 2 * Lmap->Width(), height = 2 * Lmap->Height();
        std::unique_ptr<Float[]> img(new Float[width * height]);
        float fwidth = 0.5f / std::min(width, height);
        ParallelFor(
                    [&](int64_t v) {
                        Float vp = (v + .5f) / (Float)height;
                        Float sinTheta = std::sin(Pi * (v + .5f) / height);
                        for (int u = 0; u < width; ++u) {
                            Float up = (u + .5f) / (Float)width;
                            img[u + v * width] = Lmap->Lookup(Point2f(up, vp), fwidth).y();
                            img[u + v * width] *= sinTheta;
                        }
                    },
                    height, 32);
        
        // Compute sampling distributions for rows and columns of image
        distribution.reset(new Distribution2D(img.get(), width, height));
        
        if (!Inside(this->worldSpacePosition, this->proxyVolume)) {
            Error("Light probe not inside proxy volume");
        }
        
    }
    
    Float LightProbe::ComputeInfluenceWeight(const Point3f& point) const {
        Float pointRadius = this->influenceArea.evaluateRadius(point - this->worldSpacePosition);
        
        Float clampedRadius = Clamp(pointRadius, this->influenceArea.innerRadius, this->influenceArea.outerRadius);
        Float weight = 1.0 - (clampedRadius - this->influenceArea.innerRadius) / (this->influenceArea.outerRadius - this->influenceArea.innerRadius);
        
        return weight;
    }


    Spectrum LightProbe::Power(Float worldRadius) const {
        return Pi * worldRadius * worldRadius *
        Spectrum(Lmap->Lookup(Point2f(.5f, .5f), .5f),
                 SpectrumType::Illuminant);
    }
    
    Spectrum LightProbe::Le(const Ray &ray) const {
        Float t0, t1;
        if (!this->proxyVolume.IntersectP(ray, &t0, &t1)) {
            Error("Proxy volume intersection failed for light probe.");
        }
        Float t = t0 > 0 ? t0 : t1;
        Point3f proxyVolumeHit = ray(t);
        Vector3f parallaxCorrectedSampleRayWorldSpace = proxyVolumeHit - this->worldSpacePosition;
        
        Vector3f w = Normalize(WorldToLight(parallaxCorrectedSampleRayWorldSpace));
        Point2f st(SphericalPhi(w) * Inv2Pi, SphericalTheta(w) * InvPi);
        return Spectrum(Lmap->Lookup(st), SpectrumType::Illuminant);
    }
    
    Float LightProbe::Pdf_Li(const Interaction &, const Vector3f &w) const {
        ProfilePhase _(Prof::LightPdf);
        Vector3f wi = WorldToLight(w);
        Float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
        Float sinTheta = std::sin(theta);
        if (sinTheta == 0) return 0;
        return distribution->Pdf(Point2f(phi * Inv2Pi, theta * InvPi)) /
        (2 * Pi * Pi * sinTheta);
    }

    
    std::shared_ptr<LightProbe> CreateLightProbe(const Transform &light2world, const ParamSet &paramSet) {
        Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
        Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
        std::string texmap = paramSet.FindOneFilename("mapname", "");
        
        Vector3f ellipsoidABC =  paramSet.FindOneVector3f("ellipsoiddabc", Vector3f(1, 1, 1));
        Float innerRadius = paramSet.FindOneFloat("innerradius", 1);
        Float outerRadius = paramSet.FindOneFloat("outerradius", 1);

        LightProbe::InfluenceArea influenceArea(ellipsoidABC, innerRadius, outerRadius);
        
        Point3f minProxyBound = paramSet.FindOnePoint3f("minbound", Point3f(0, 0, 0));
        Point3f maxProxyBound = paramSet.FindOnePoint3f("maxbound", Point3f(1, 1, 1));
        
        Bounds3f proxyVolume = Bounds3f(minProxyBound, maxProxyBound);

        return std::make_shared<LightProbe>(light2world, L * sc, texmap, proxyVolume, influenceArea);
    }
    
}  // namespace pbrt
