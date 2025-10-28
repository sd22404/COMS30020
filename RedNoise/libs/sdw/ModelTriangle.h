#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "Material.h"
#include "TexturePoint.h"
#include "Vertex.h"

struct ModelTriangle {
	std::array<Vertex, 3> vertices{};
	Material material{};
	glm::vec3 normal{};
	std::string texture{};
	std::string normalMap{};

	Vertex &v0();
	Vertex &v1();
	Vertex &v2();

	ModelTriangle();
	ModelTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2, Material &mat);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
