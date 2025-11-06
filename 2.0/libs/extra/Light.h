#pragma once

#include <glm/glm.hpp>

class Light {
public:
    Light(const glm::vec3 &position, const float &intensity)
        : position(position), colour(glm::vec3(1, 1, 1)), intensity(intensity) {}
    Light(const glm::vec3 &position, glm::vec3 &colour, const float &intensity)
        : position(position), colour(colour), intensity(intensity) {}

    glm::vec3 position;
    glm::vec3 colour;
    float intensity;
};