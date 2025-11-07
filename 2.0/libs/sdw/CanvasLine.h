#pragma once
#include <vector>
#include "CanvasPoint.h"
#include "TexturePoint.h"

struct CanvasLine {
    std::vector<CanvasPoint> points{};

    bool isOffScreen(float width, float height) const;
    static std::vector<TexturePoint> interpolateTexturePoints(const TexturePoint &from, const TexturePoint &to, long numberOfValues);
    CanvasLine(const CanvasPoint &from, const CanvasPoint &to, long numberOfValues);
};
