#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, Colour trigColour) :
		vertices({{v0, v1, v2}}), texturePoints(), colour(std::move(trigColour)), normal() {}

glm::vec3 &ModelTriangle::v0() {
	return {vertices[0].position};
}

glm::vec3 &ModelTriangle::v1() {
	return {vertices[1].position};
}

glm::vec3 &ModelTriangle::v2() {
	return {vertices[2].position};
}


std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
	os << "(" << triangle.vertices[0] << ")\n";
	os << "(" << triangle.vertices[1] << ")\n";
	os << "(" << triangle.vertices[2] << ")\n";
	return os;
}
