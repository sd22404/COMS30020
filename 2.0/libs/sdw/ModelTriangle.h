#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include <TexturePoint.h>
#include "Vertex.h"
#include "Material.h"

struct ModelTriangle {
	std::array<Vertex, 3> vertices{};
	glm::vec3 normal{};

	Vertex operator[](size_t i) const;
	Vertex &operator[](size_t i);

	ModelTriangle();
	ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
