#pragma once

#include <glm/glm.hpp>
#include <string>
#include <utility>
#include "TextureMap.h"

struct Material {
    std::string name;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
    float ambient;
    float reflectivity;

    TextureMap texture = TextureMap();
    TextureMap normalMap = TextureMap();

    // bool emissive;
    bool glassy;
    float refractiveIndex;

    Material() : name("default"), diffuse(glm::vec3(1, 1, 1)), specular(1, 1, 1), shininess(256.0f), ambient(0.1f), reflectivity(0) {}
    Material(std::string name, const glm::vec3 colour) : name(std::move(name)), diffuse(colour), specular(1, 1, 1), shininess(256.0f), ambient(0.1f), reflectivity(0) {}
};