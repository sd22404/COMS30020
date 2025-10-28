#pragma once

#include <iostream>
#include <glm/detail/type_vec3.hpp>

struct Material {
	std::string name;
	glm::vec3 diffuse;
	bool emissive{};
	bool mirrored{};
	bool glassy{};
	float refractiveIndex{};

	Material();
	explicit Material(glm::vec3 colour);
	Material(std::string n, glm::vec3 colour);
};

std::ostream &operator<<(std::ostream &os, const Material &colour);
