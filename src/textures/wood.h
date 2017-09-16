
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

#ifndef PBRT_TEXTURES_WOOD_H
#define PBRT_TEXTURES_WOOD_H

// textures/wood.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// Original inspiration: https://www.shadertoy.com/view/ldscDM

namespace pbrt {
    inline Float repramp(Float x) {
        return pow(sin(x)*0.5+0.5, 8.0) + cos(x)*0.7 + 0.7;
    }

// WoodTexture Declarations
class WoodTexture : public Texture<Spectrum> {
  public:
    // WoodTexture Public Methods
    WoodTexture(std::unique_ptr<TextureMapping3D> mapping, const Spectrum& woodColour, const Spectrum& ringColour, const Float ringNoisiness, const Float ringDetail, const Float ringSpacing, const Float fleckIntensity)
    : mapping(std::move(mapping)), woodColour(woodColour), ringColour(ringColour), ringNoisiness(ringNoisiness), ringDetail(ringDetail), ringSpacing(ringSpacing), fleckIntensity(fleckIntensity) {}
    
    Spectrum Evaluate(const SurfaceInteraction &si) const {
        Vector3f dpdx, dpdy;
        Point3f p = mapping->Map(si, &dpdx, &dpdy);
        
        Point3f scaledPoint(p.x * 20, p.y * 20, p.z * 20);
        
        auto noise1 = Noise(p.x * 8, p.y * 1.5, p.z * 8);
        auto noise2 = Noise(-p.x * 8 + 4.5678, -p.y * 1.5 + 4.5678, -p.z * 8 + 4.5678);
        auto noise3 = Noise(p.x * 64.0, p.y * 12.8, p.z * 64);
        
        auto posYX = Vector2f(p.y, p.x);
        
        Float rings = repramp((posYX + Vector2f(noise1, noise2) * 0.1 * ringNoisiness).Length() * 128.0 * ringSpacing) / 1.8;
        rings -= Noise(p.x, p.y, p.z) * 0.75;
        
        Spectrum texColour = Lerp(rings, ringColour, woodColour);
        texColour = texColour.Clamp();
        
        // Add a small amount of noise to avoid a uniform colour.
        float noisyRoughness = noise3 * 0.1 + 0.9;
        texColour *= noisyRoughness;
        
        // Add detail noise to the rings:
        Float detail = Turbulence(scaledPoint, dpdx, dpdy, 0.9, 12);
        texColour = texColour + 0.1 * ringDetail * (1 - rings) * Spectrum(detail - 1.75);
        
        // Add dark lines
        Point3f linesP(p.x * 0.08 - 0.2 * rings, p.y - noise1, p.z + 0.2 * noise2);
        Spectrum fleckColour = texColour;
        if (FBm(linesP, dpdx, dpdy, 1.5, 11) > 0.995) {
            fleckColour = Lerp(1 - (rings * 0.7 + 0.3), ringColour, woodColour) * 0.5;
        }
        texColour = Lerp(this->fleckIntensity, texColour, fleckColour);
        
        texColour = texColour.Clamp(0, 1);
        
        return texColour;
    }

  private:
    // WoodTexture Private Data
    std::unique_ptr<TextureMapping3D> mapping;
    const Spectrum woodColour, ringColour;
    const Float ringNoisiness, ringDetail, ringSpacing, fleckIntensity;
};

WoodTexture *CreateWoodSpectrumTexture(
    const Transform &tex2world, const TextureParams &tp);

}  // namespace pbrt

#endif  // PBRT_TEXTURES_WOOD_H
