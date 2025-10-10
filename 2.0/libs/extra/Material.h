#pragma once

#include <glm/glm.hpp>
#include <string>

class Material {
public:
    std::string name;
    glm::vec3 colour;
    bool mirrored;
    bool emissive;
    bool glassy;
    float refractiveIndex;

    Material() : name("default"), colour(glm::vec3(1, 1, 1)), mirrored(false), emissive(false), glassy(false), refractiveIndex(1.0f) {}
    Material(std::string name, glm::vec3 colour) : name(name), colour(colour), mirrored(false), emissive(false), glassy(false), refractiveIndex(1.0f) {}
};