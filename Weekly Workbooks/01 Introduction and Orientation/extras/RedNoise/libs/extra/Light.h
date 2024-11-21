#pragma once

#include <ostream>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

#define DEFAULT_INT 25.0f
#define DEFAULT_POS glm::vec3(0.0f, 0.9f, 0.0f)

struct Light {
    glm::vec3 defaultPosition{DEFAULT_POS};
    glm::vec3 position{defaultPosition};
    float intensity{DEFAULT_INT};

    Light();
    explicit Light(float intensity);
    Light(float intensity, glm::vec3 position);

    friend std::ostream &operator<<(std::ostream &os, const Light &light);
};