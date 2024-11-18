#pragma once

#include <ostream>
#include <glm/detail/type_mat3x3.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

struct Camera {
    glm::vec3 defaultPosition{};
    glm::mat3 defaultRotation{};
    float defaultFocal{};
    glm::vec3 position{};
    glm::mat3 rotation{};
    float focalLength{};
    float moveSpeed{};
    float altMoveSpeed{};
    float lookSpeed{};
    float altLookSpeed{};

    Camera();
    Camera(float focalLength);
    Camera(glm::vec3 position, glm::mat3 rotation);
    Camera(glm::vec3 position, glm::mat3 rotation, float focalLength);

    friend std::ostream &operator<<(std::ostream &os, const Camera &cam);
};