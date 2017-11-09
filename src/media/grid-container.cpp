
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


// media/grid.cpp*
#include "media/grid-container.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"
#include "interaction.h"

namespace pbrt {

STAT_RATIO("Media/Grid steps per Tr() call", nTrSteps, nTrCalls);

Spectrum GridDensityMediaContainer::Sample(const Ray &rWorld, Sampler &sampler,
                                   MemoryArena &arena,
                                   MediumInteraction *mi) const {
    const Float firstGridTime = 0.0; // hardcoded to 0 for now.
    Float time = rWorld.time;
    
    if (this->gridDensityMedia.size() == 1 || time <= firstGridTime) {
        return this->gridDensityMedia.front()->Sample(rWorld, sampler, arena, mi);
    }
    
    if (time >= this->endTime) {
        return this->gridDensityMedia.back()->Sample(rWorld, sampler, arena, mi);
    }
    
    
    int endGridIndex = this->gridDensityMedia.size() - 1;
    for (size_t i = 0; i < this->gridDensityMedia.size(); i++) {
        Float frameTime = i * frameInterval;
        
        if (frameTime > time) {
            endGridIndex = i;
            break;
        }
    }
    
    int startGridIndex = endGridIndex == 0 ? endGridIndex : endGridIndex - 1;
    
    Spectrum startGridSpectrum = this->gridDensityMedia[startGridIndex]->Sample(rWorld, sampler, arena, mi);
    if (startGridIndex == endGridIndex) {
        return startGridSpectrum;
    }
    
    Float startGridTime = startGridIndex * frameInterval;
    Float endGridTime = endGridIndex * frameInterval;
    
    Float dt = (time - startGridTime) / (endGridTime - startGridTime);
    Spectrum endGridSpectrum = this->gridDensityMedia[endGridIndex]->Sample(rWorld, sampler, arena, mi);
    
    // Interpolate translation at _dt_
    return (1 - dt) * startGridSpectrum + dt * endGridSpectrum;
}

Spectrum GridDensityMediaContainer::Tr(const Ray &rWorld, Sampler &sampler) const {
    const Float firstGridTime = 0.0; // hardcoded to 0 for now.
    Float time = rWorld.time;
    
    if (this->gridDensityMedia.size() == 1 || time <= firstGridTime) {
        return this->gridDensityMedia.front()->Tr(rWorld, sampler);
    }
    
    if (time >= this->endTime) {
        return this->gridDensityMedia.back()->Tr(rWorld, sampler);
    }
    
    
    int endGridIndex = this->gridDensityMedia.size() - 1;
    for (size_t i = 0; i < this->gridDensityMedia.size(); i++) {
        Float frameTime = i * frameInterval;
        
        if (frameTime > time) {
            endGridIndex = i;
            break;
        }
    }
    
    int startGridIndex = endGridIndex == 1 ? endGridIndex : endGridIndex - 1;
    
    Spectrum startGridSpectrum = this->gridDensityMedia[startGridIndex]->Tr(rWorld, sampler);
    if (startGridIndex == endGridIndex) {
        return startGridSpectrum;
    }
    
    Float startGridTime = startGridIndex * frameInterval;
    Float endGridTime = endGridIndex * frameInterval;
    
    Float dt = (time - startGridTime) / (endGridTime - startGridTime);
    
    Spectrum endGridSpectrum = this->gridDensityMedia[endGridIndex]->Tr(rWorld, sampler);
    
    // Interpolate translation at _dt_
    return (1 - dt) * startGridSpectrum + dt * endGridSpectrum;
}

}  // namespace pbrt
