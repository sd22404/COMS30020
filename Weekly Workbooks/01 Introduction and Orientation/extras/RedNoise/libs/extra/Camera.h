#pragma once

#include <ostream>
#include <glm/detail/type_mat3x3.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

#define AMBIENT 0.1f

enum Render { WIRE, RASTER, RAYTRACE };
enum Shading { FLAT, GOURAUD, PHONG };

struct Camera {
    glm::vec3 defaultPosition{glm::vec3(0.0f, 0.0f, 4.0f)};
    glm::mat3 defaultRotation{glm::mat3(glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f))};
    glm::vec3 position{defaultPosition};
    glm::mat3 rotation{defaultRotation};
    float focalLength{2.0f};
    float moveSpeed{0.05f};
    float altMoveSpeed{0.25f};
    float lookSpeed{0.5f};
    float altLookSpeed{2.5f};

    Render renderMode{RASTER};
    Shading shadingMode{FLAT};
    bool orbit{};
    float ambientLight{AMBIENT};

    Camera();
    explicit Camera(float focalLength);
    Camera(glm::vec3 position, float focalLength);
    Camera(glm::vec3 position, glm::mat3 rotation, float focalLength);

    friend std::ostream &operator<<(std::ostream &os, const Camera &cam);
};