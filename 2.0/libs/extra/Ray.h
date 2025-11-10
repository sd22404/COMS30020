#pragma once

#include "glm/vec3.hpp"

struct Ray {
    glm::vec3 start;
    glm::vec3 dir;
    float dist;
    Ray() = default;
    Ray(const glm::vec3 start, const glm::vec3 direction) : start(start), dir(direction), dist(0) {};
    Ray(const glm::vec3 start, const glm::vec3 direction, const float distance) : start(start), dir(direction), dist(distance) {};
};
