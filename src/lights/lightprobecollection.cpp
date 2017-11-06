
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
#include "lights/lightprobecollection.h"
#include "lights/lightprobe.h"
#include "imageio.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {
    
    // LightProbe Method Definitions
    LightProbeCollection::LightProbeCollection(const std::vector<std::shared_ptr<LightProbe>> lightProbes)
    : Light((int)LightFlags::Infinite, Transform(), MediumInterface(), 1), lightProbes(lightProbes) {
        
    
    }
    
    Spectrum LightProbeCollection::Power() const {
        Spectrum powerSum(0.0);
        for (auto& lightProbe : this->lightProbes) {
            powerSum += lightProbe->Power(this->worldRadius);
        }
        return powerSum/this->lightProbes.size();
    }
    
    Spectrum LightProbeCollection::Le(const RayDifferential &ray) const {
        Float weights[this->lightProbes.size()];
        int maxWeightIndex;
        Float weightSum = this->ComputeLightProbeBlendWeights(&weights[0], ray.o, maxWeightIndex);
        
        Spectrum Le(0.0);
        Float invWeightSum = 1.0/weightSum;
        for (size_t i = 0; i < this->lightProbes.size(); i++) {
            if (weights[i] != 0) {
                Le += weights[i] * invWeightSum * this->lightProbes[i]->Le(ray);
            }
        }

        return Le;
    }
    
    Spectrum LightProbeCollection::Sample_Li(const Interaction &ref, const Point2f &u,
                                   Vector3f *wi, Float *pdf,
                                   VisibilityTester *vis) const {
        
        Float weights[this->lightProbes.size()];
        int maxWeightIndex;
        Float weightSum = this->ComputeLightProbeBlendWeights(&weights[0], ref.p, maxWeightIndex);
        
        const std::shared_ptr<LightProbe>& maxWeightProbe = this->lightProbes[maxWeightIndex];
        
        // Find $(u,v)$ sample coordinates in infinite light texture
        Float mapPdf;
        Point2f uv = maxWeightProbe->distribution->SampleContinuous(u, &mapPdf);
        if (mapPdf == 0) return Spectrum(0.f);
        
        // Convert infinite light sample point to direction
        Float theta = uv[1] * Pi, phi = uv[0] * 2 * Pi;
        Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
        Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
        *wi = maxWeightProbe->LightToWorld(Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta));
        
        // Compute PDF for sampled infinite light direction
        *pdf = mapPdf / (2 * Pi * Pi * sinTheta);
        if (sinTheta == 0) *pdf = 0;
        
        // Return radiance value for infinite light direction
        *vis = VisibilityTester(ref, Interaction(ref.p + *wi * (2 * worldRadius),
                                                 ref.time, mediumInterface));
        
        Ray ray(ref.p, *wi);
        
        Float invWeightSum = 1.0/weightSum;
        
        Spectrum Le(0.0);
        for (size_t i = 0; i < this->lightProbes.size(); i++) {
            if (weights[i] != 0) {
                Le += weights[i] * invWeightSum * this->lightProbes[i]->Le(ray);
            }
        }
        
        return Le;
    }
    
    Float LightProbeCollection::Pdf_Li(const Interaction &i, const Vector3f &w) const {
        ProfilePhase _(Prof::LightPdf);
        
        Float weights[this->lightProbes.size()];
        int maxWeightIndex;
        this->ComputeLightProbeBlendWeights(&weights[0], i.p, maxWeightIndex);
        
        return this->lightProbes[maxWeightIndex]->Pdf_Li(i, w);
    }
    
    Spectrum LightProbeCollection::Sample_Le(const Point2f &u1, const Point2f &u2,
                                   Float time, Ray *ray, Normal3f *nLight,
                                   Float *pdfPos, Float *pdfDir) const {
        Error("Not Supported");
        return 0;
    }
    
    void LightProbeCollection::Pdf_Le(const Ray &ray, const Normal3f &, Float *pdfPos,
                            Float *pdfDir) const {
        Error("Not Supported");
        return;
    }
    

    
    Float LightProbeCollection::ComputeLightProbeBlendWeights(Float *weights, const Point3f& intersection, int& maxWeightIndex) const {
        for (size_t i = 0; i < this->lightProbes.size(); i++) {
            weights[i] = this->lightProbes[i]->ComputeInfluenceWeight(intersection);
        }
        
        Float totalWeight = weights[0];
        maxWeightIndex = 0;
        for (size_t i = 1; i < this->lightProbes.size(); i++) {
            Float weight = weights[i];
            totalWeight += weight;
            if (weight > weights[maxWeightIndex]) {
                maxWeightIndex = i;
            }
        }
    
        return totalWeight;
    }
    
    std::shared_ptr<LightProbeCollection> CreateLightProbeCollection(const std::vector<std::shared_ptr<LightProbe>> lightProbes) {
        return std::make_shared<LightProbeCollection>(lightProbes);
    }
    
}  // namespace pbrt

