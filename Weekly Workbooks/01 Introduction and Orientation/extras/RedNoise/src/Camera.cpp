#include "Camera.h"

#define DEFAULT_POS glm::vec3(0.0f, 0.0f, 4.0f)
#define DEFAULT_ROT glm::mat3(glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f))
#define DEFAULT_FOCAL 2.0f
#define MOVE_SPEED 0.05f

Camera::Camera() :
    defaultPosition(DEFAULT_POS),
    defaultRotation(DEFAULT_ROT),
    defaultFocal(DEFAULT_FOCAL),
    position(defaultPosition),
    rotation(defaultRotation),
    focalLength(defaultFocal),
    moveSpeed(MOVE_SPEED) {}

Camera::Camera(float focalLength) :
    defaultPosition(DEFAULT_POS),
    defaultRotation(DEFAULT_ROT),
    defaultFocal(focalLength),
    position(defaultPosition),
    rotation(defaultRotation),
    focalLength(defaultFocal),
    moveSpeed(MOVE_SPEED) {}

Camera::Camera(glm::vec3 position, glm::mat3 rotation) :
    defaultPosition(position),
    defaultRotation(rotation),
    defaultFocal(DEFAULT_FOCAL),
    position(defaultPosition),
    rotation(defaultRotation),
    focalLength(defaultFocal),
    moveSpeed(MOVE_SPEED) {}

Camera::Camera(glm::vec3 position, glm::mat3 rotation, float focalLength) :
    defaultPosition(position),
    defaultRotation(rotation),
    defaultFocal(focalLength),
    position(defaultPosition),
    rotation(defaultRotation),
    focalLength(defaultFocal),
    moveSpeed(MOVE_SPEED) {}

std::ostream &operator<<(std::ostream &os, const Camera &cam) {
    os << "(" << to_string(cam.position) << ", " << to_string(cam.rotation) << ") " << cam.focalLength;
    return os;
}