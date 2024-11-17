#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, Colour trigColour) :
		vertices({{v0, v1, v2}}), colour(std::move(trigColour)), normal(), texture() {}

Vertex &ModelTriangle::v0() {
	return {vertices[0]};
}

Vertex &ModelTriangle::v1() {
	return {vertices[1]};
}

Vertex &ModelTriangle::v2() {
	return {vertices[2]};
}


std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
	os << "(" << triangle.vertices[0] << ")\n";
	os << "(" << triangle.vertices[1] << ")\n";
	os << "(" << triangle.vertices[2] << ")\n";
	return os;
}
