#include "CanvasTriangle.h"

bool CanvasTriangle::isOffScreen(const float width, const float height) const {
    return (vertices[0].isOffScreen(width, height) &&
            vertices[1].isOffScreen(width, height) &&
            vertices[2].isOffScreen(width, height));
}

CanvasTriangle::CanvasTriangle() = default;
CanvasTriangle::CanvasTriangle(const CanvasPoint &v0, const CanvasPoint &v1, const CanvasPoint &v2) :
    vertices({{v0, v1, v2}}) {}

CanvasPoint CanvasTriangle::operator[](const size_t i) const {
    return vertices[i];
}

CanvasPoint &CanvasTriangle::operator[](const size_t i) {
    return vertices[i];
}

std::ostream &operator<<(std::ostream &os, const CanvasTriangle &triangle) {
	os << triangle[0] << triangle[1] << triangle[2];
	return os;
}
