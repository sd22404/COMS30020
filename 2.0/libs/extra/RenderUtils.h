#pragma once

#include <vector>
#include "CanvasPoint.h"
#include "TexturePoint.h"

class RenderUtils {
public:
    static std::vector<CanvasPoint> interpolateCanvasPoints(CanvasPoint from, CanvasPoint to, float numberOfValues);
    static std::vector<TexturePoint> interpolateTexturePoints(TexturePoint from, TexturePoint to, float numberOfValues);
};