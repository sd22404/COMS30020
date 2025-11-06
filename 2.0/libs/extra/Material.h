#pragma once

#include <glm/glm.hpp>
#include <string>

class Material {
public:
    std::string name;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
    float ambient;
    float reflectivity;
    bool emissive;
    bool glassy;
    float refractiveIndex;

    Material() : name("default"), diffuse(glm::vec3(1, 1, 1)), specular(1, 1, 1), shininess(32.0f), ambient(0.1f), reflectivity(0), emissive(false), glassy(false), refractiveIndex(1.0f) {}
    Material(std::string name, glm::vec3 colour) : name(name), diffuse(colour), specular(1, 1, 1), shininess(32.0f), ambient(0.1f), reflectivity(0), emissive(false), glassy(false), refractiveIndex(1.0f) {}
};