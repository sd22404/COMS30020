#pragma once

#include "glm/vec3.hpp"

class Ray {
public:
    glm::vec3 start;
    glm::vec3 dir;
    float dist;
    Ray(glm::vec3 start, glm::vec3 direction) : start(start), dir(direction) {};
    Ray(glm::vec3 start, glm::vec3 direction, float distance) : start(start), dir(direction), dist(distance) {};
};
