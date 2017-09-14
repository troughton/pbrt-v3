
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


// textures/wood.cpp*
#include "textures/wood.h"

namespace pbrt {
    
WoodTexture*CreateWoodSpectrumTexture(
    const Transform &tex2world, const TextureParams &tp) {
    // Initialize 3D texture mapping _map_ from _tp_
    std::unique_ptr<TextureMapping3D> map(new IdentityMapping3D(tex2world));
    
    Float colour1RGB[3] = { 0.43, 0.27, 0.107 };
    Float colour2RGB[3] = { 0.6, 0.44, 0.196 };
    Spectrum colourA = Spectrum::FromRGB(colour1RGB);
    Spectrum colourB = Spectrum::FromRGB(colour2RGB);
    
    return new WoodTexture(std::move(map),
                                       tp.FindSpectrum("colourA", colourA),
                           tp.FindSpectrum("colourB", colourB));
}

}  // namespace pbrt
