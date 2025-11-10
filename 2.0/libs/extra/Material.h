#pragma once

#include <glm/glm.hpp>
#include <string>
#include <utility>
#include "TextureMap.h"

struct Material {
    std::string name = "default";
    glm::vec3 diffuse = glm::vec3(1);
    glm::vec3 specular = glm::vec3(1);
    float shininess = 256.0f;
    float ambient = 0.1f;
    float reflectivity = 0.0f;
    float transparency = 0.0f;
    float refractiveIndex = 1.0f;
    bool emissive = false;
    TextureMap texture = TextureMap();
    TextureMap normalMap = TextureMap();

    Material() = default;
    Material(std::string name, const glm::vec3 colour) : name(std::move(name)), diffuse(colour) {}
};