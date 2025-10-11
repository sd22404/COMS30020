#include "CanvasLine.h"

std::vector<TexturePoint> CanvasLine::interpolateTexturePoints(TexturePoint &from, TexturePoint &to, float numberOfValues) {
    // if one or fewer values, return only the start point (as it will be the same as the end point)
    std::vector<TexturePoint> result = {from};
    if (numberOfValues < 1) return result;

    TexturePoint interval = {to.x - from.x, to.y - from.y};
    float stepX = interval.x / numberOfValues;
    float stepY = interval.y / numberOfValues;

    for (float i = 1; i < numberOfValues; i++) {
        TexturePoint point = {from.x + i * stepX, from.y + i * stepY};
        result.push_back(point);
    }
    return result;
}

bool CanvasLine::isOffScreen(float width, float height) {
    return (points.front().x >= width && points.back().x >= width) ||
    (points.front().y >= height && points.back().y >= height) ||
    (points.front().x < 0 && points.back().x < 0) ||
    (points.front().y < 0 && points.back().y < 0) ||
    (points.front().depth < 0 && points.back().depth < 0);
}

CanvasLine::CanvasLine(CanvasPoint &from, CanvasPoint &to, float numberOfValues) {
    // if one or fewer values, return only the start point (as it will be the same as the end point)
    points = {from};
    if (numberOfValues < 1) return;
    CanvasPoint interval = {to.x - from.x, to.y - from.y, to.depth - from.depth};
    float stepX = interval.x / numberOfValues;
    float stepY = interval.y / numberOfValues;
    float stepD = interval.depth / numberOfValues;

    std::vector<TexturePoint> tps = interpolateTexturePoints(from.texturePoint, to.texturePoint, numberOfValues);
    for (float i = 1; i < numberOfValues; i++) {
        CanvasPoint point = {from.x + i * stepX, from.y + i * stepY, from.depth + i * stepD};
        point.texturePoint = tps[i];
        points.push_back(point);
    }
}