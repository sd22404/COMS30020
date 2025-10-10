#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "TexturePoint.h"
#include "../extra/Vertex.h"
#include "../extra/Material.h"

struct ModelTriangle {
	std::array<Vertex, 3> vertices{};
	std::array<TexturePoint, 3> texturePoints{};
	Material material{};
	glm::vec3 normal{};
	std::string texture;
	std::string normalMap;

	Vertex &v0() { return vertices[0]; }
	Vertex &v1() { return vertices[1]; }
	Vertex &v2() { return vertices[2]; }

	ModelTriangle();
	ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, Material &mat);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
