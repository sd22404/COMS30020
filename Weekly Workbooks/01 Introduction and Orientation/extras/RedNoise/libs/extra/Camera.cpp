#include "Camera.h"



Camera::Camera() = default;

Camera::Camera(float focalLength) : focalLength(focalLength) {}

Camera::Camera(glm::vec3 position, glm::mat3 rotation) :
    defaultPosition(position),
    defaultRotation(rotation),
    position(defaultPosition),
    rotation(defaultRotation) {}

Camera::Camera(glm::vec3 position, glm::mat3 rotation, float focalLength) :
    defaultPosition(position),
    defaultRotation(rotation),
    position(defaultPosition),
    rotation(defaultRotation),
    focalLength(focalLength) {}

std::ostream &operator<<(std::ostream &os, const Camera &cam) {
    os << "(" << to_string(cam.position) << ", " << to_string(cam.rotation) << ") " << cam.focalLength;
    return os;
}