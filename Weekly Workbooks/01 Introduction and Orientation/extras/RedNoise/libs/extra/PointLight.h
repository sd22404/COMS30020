#pragma once

#include <ostream>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>


struct PointLight {
    glm::vec3 defaultPosition{};
    glm::vec3 position{};
    float intensity{};

    PointLight();
    PointLight(float intensity);
    PointLight(float intensity, glm::vec3 position);

    friend std::ostream &operator<<(std::ostream &os, const PointLight &light);
};