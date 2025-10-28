#include "Material.h"
#include <utility>

Material::Material() = default;
Material::Material(glm::vec3 col) : diffuse(col) {}
Material::Material(std::string n, glm::vec3 col) : name(std::move(n)), diffuse(col) {};

std::ostream &operator<<(std::ostream &os, const Material &colour) {
	os << colour.name << " ["
	   << colour.diffuse.r << ", "
	   << colour.diffuse.g << ", "
	   << colour.diffuse.b << "]";
	return os;
}