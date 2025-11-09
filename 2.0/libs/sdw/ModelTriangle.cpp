#include "ModelTriangle.h"

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2) :
		vertices({{v0, v1, v2}}), normal() {}

Vertex ModelTriangle::operator[](size_t i) const {
	return vertices[i];
}

Vertex &ModelTriangle::operator[](size_t i) {
	return vertices[i];
}

std::ostream &operator<<(std::ostream &os, ModelTriangle &triangle) {
	os << "(" << triangle[0].position.x << ", " << triangle[0].position.y << ", " << triangle[0].position.z << ")\n";
	os << "(" << triangle[1].position.x << ", " << triangle[1].position.y << ", " << triangle[1].position.z << ")\n";
	os << "(" << triangle[2].position.x << ", " << triangle[2].position.y << ", " << triangle[2].position.z << ")\n";
	return os;
}
