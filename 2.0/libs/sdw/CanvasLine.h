#pragma once
#include <vector>
#include "CanvasPoint.h"
#include "TexturePoint.h"

struct CanvasLine {
    std::vector<CanvasPoint> points{};

    bool isOffScreen(float width, float height);
    static std::vector<TexturePoint> interpolateTexturePoints(TexturePoint &from, TexturePoint &to, float numberOfValues);
    CanvasLine(CanvasPoint &from, CanvasPoint &to, float numberOfValues);
};
