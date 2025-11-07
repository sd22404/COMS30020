#include "CanvasLine.h"

std::vector<TexturePoint> CanvasLine::interpolateTexturePoints(const TexturePoint &from, const TexturePoint &to, const long numberOfValues) {
    // if one or fewer values, return only the start point (as it will be the same as the end point)
    std::vector<TexturePoint> result = {from};
    if (numberOfValues < 1) return result;

    TexturePoint interval = {to.x - from.x, to.y - from.y};
    const float stepX = interval.x / static_cast<float>(numberOfValues);
    const float stepY = interval.y / static_cast<float>(numberOfValues);

    for (int i = 1; i < numberOfValues; i++) {
        TexturePoint point = {from.x + static_cast<float>(i) * stepX, from.y + static_cast<float>(i) * stepY};
        result.push_back(point);
    }
    return result;
}

bool CanvasLine::isOffScreen(const float width, const float height) const {
    return (points.front().x >= width && points.back().x >= width) ||
    (points.front().y >= height && points.back().y >= height) ||
    (points.front().x < 0 && points.back().x < 0) ||
    (points.front().y < 0 && points.back().y < 0) ||
    (points.front().depth < 0 && points.back().depth < 0);
}

CanvasLine::CanvasLine(const CanvasPoint &from, const CanvasPoint &to, const long numberOfValues) {
    // if one or fewer values, return only the start point (as it will be the same as the end point)
    points = {from};
    if (numberOfValues < 1) return;
    const auto interval = CanvasPoint(to.x - from.x, to.y - from.y, to.depth - from.depth);
    const float stepX = interval.x / static_cast<float>(numberOfValues);
    const float stepY = interval.y / static_cast<float>(numberOfValues);
    const float stepD = interval.depth / static_cast<float>(numberOfValues);

    const std::vector<TexturePoint> tps = interpolateTexturePoints(from.texturePoint, to.texturePoint, numberOfValues);
    for (int i = 1; i < numberOfValues; i++) {
        CanvasPoint point = {from.x + static_cast<float>(i) * stepX, from.y + static_cast<float>(i) * stepY, from.depth + static_cast<float>(i) * stepD};
        point.texturePoint = tps[i];
        points.push_back(point);
    }
}