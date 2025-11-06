#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "TexturePoint.h"
#include "../extra/Vertex.h"
#include "../extra/Material.h"

struct ModelTriangle {
	std::array<Vertex, 3> vertices{};
	Material material{};
	glm::vec3 normal{};
	std::string texture;
	std::string normalMap;

	Vertex operator[](size_t i) const;
	Vertex &operator[](size_t i);

	ModelTriangle();
	ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, Material &mat);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
