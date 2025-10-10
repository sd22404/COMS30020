#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, Material &mat) :
		vertices({{v0, v1, v2}}), texturePoints(), material(mat), normal() {}

std::ostream &operator<<(std::ostream &os, ModelTriangle &triangle) {
	os << "(" << triangle.v0().position.x << ", " << triangle.v0().position.y << ", " << triangle.v0().position.z << ")\n";
	os << "(" << triangle.v1().position.x << ", " << triangle.v1().position.y << ", " << triangle.v1().position.z << ")\n";
	os << "(" << triangle.v2().position.x << ", " << triangle.v2().position.y << ", " << triangle.v2().position.z << ")\n";
	return os;
}
