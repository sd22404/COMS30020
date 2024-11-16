#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "Colour.h"
#include "TexturePoint.h"
#include "../extra/Vertex.h"

struct ModelTriangle {
	std::array<Vertex, 3> vertices{};
	std::array<TexturePoint, 3> texturePoints{};
	std::string texture{};
	Colour colour{};
	glm::vec3 normal{};

	Vertex &v0();
	Vertex &v1();
	Vertex &v2();

	ModelTriangle();
	ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, Colour trigColour);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
