#pragma once

#include <ostream>
#include <glm/detail/type_mat3x3.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

#define DEFAULT_POS glm::vec3(0.0f, 0.0f, 4.0f)
#define DEFAULT_ROT glm::mat3(glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f))
#define DEFAULT_FOCAL 2.0f
#define DEFAULT_MOVE_SPEED 0.05f
#define FAST_MOVE_SPEED 0.25f
#define DEFAULT_LOOK_SPEED 0.5f
#define FAST_LOOK_SPEED 2.5f

struct Camera {
    glm::vec3 defaultPosition{DEFAULT_POS};
    glm::mat3 defaultRotation{DEFAULT_ROT};
    glm::vec3 position{defaultPosition};
    glm::mat3 rotation{defaultRotation};
    float focalLength{DEFAULT_FOCAL};
    float moveSpeed{DEFAULT_MOVE_SPEED};
    float altMoveSpeed{FAST_MOVE_SPEED};
    float lookSpeed{DEFAULT_LOOK_SPEED};
    float altLookSpeed{FAST_LOOK_SPEED};

    Camera();
    explicit Camera(float focalLength);
    Camera(glm::vec3 position, glm::mat3 rotation);
    Camera(glm::vec3 position, glm::mat3 rotation, float focalLength);

    friend std::ostream &operator<<(std::ostream &os, const Camera &cam);
};