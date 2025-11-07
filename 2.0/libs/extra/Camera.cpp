#include "Camera.h"

CanvasPoint Camera::projectVertex(const Vertex &vertex, const float canvasScale) const {
    // vertex in terms of camera coordinates
	const glm::vec3 finalPos = (vertex.position - position) * rotation;
	// transform onto image plane
	const float u = -canvasScale * focalLength * finalPos.x / finalPos.z + static_cast<float>(width) / 2.0f;
	const float v = canvasScale * focalLength * finalPos.y / finalPos.z + static_cast<float>(height) / 2.0f;
    CanvasPoint point(u, v, -1 / finalPos.z);
    point.texturePoint = vertex.texturePoint;
	return point;
}

Ray Camera::projectRay(const int x, const int y, const float canvasScale) const {
    // convert from SDL coordinate system into 3D/model coordinate system
    const float sdlX = static_cast<float>(x) - static_cast<float>(width) / 2.0f;
    const float sdlY = static_cast<float>(-y) + static_cast<float>(height) / 2.0f;
    // generate canvas point in 3D space, adjusted by cameraOrientation
    const glm::vec3 canvasPoint3D = position + glm::vec3(sdlX, sdlY, -(focalLength * canvasScale)) * inverse(rotation);
    // subtract cameraPosition and normalise to get ray direction
    return {position, normalize(canvasPoint3D - position)};
}

glm::mat3 Camera::rotateY(const float angle) {
    return {
        glm::vec3(std::cos(angle), 0, -std::sin(angle)),
        glm::vec3(0, 1, 0),
        glm::vec3(std::sin(angle), 0, std::cos(angle))
    };
}

float Camera::degToRad(const float deg) {
    return M_PI * deg / 180;
}

void Camera::lookAt(const glm::vec3 target) {
    const glm::vec3 forward = normalize(position - target);
    const glm::vec3 right = normalize(cross(glm::vec3(0, 1, 0), forward));
    const glm::vec3 up = normalize(cross(forward, right));
    rotation = glm::mat3(right, up, forward);
}

void Camera::move(const Direction dir) {
    switch (dir) {
        case UP:
            position.y += speed;
            break;
        case DOWN:
            position.y -= speed;
            break;
        case LEFT:
            position.x -= speed;
            break;
        case RIGHT:
            position.x += speed;
            break;
        case FORWARD:
            position.z -= speed;
            break;
        case BACKWARD:
            position.z += speed;
            break;
        default:
            break;
    }
}

void Camera::reset() {
    position = startPosition;
    rotation = glm::mat3(glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec3(0, 0, 1));
}

void Camera::toggleOrbit() { orbiting = !orbiting; }

void Camera::orbit() {
    if (!orbiting) return;
    position = position * rotateY(degToRad(speed * 10.0f));
    lookAt({0, 0, 0});
}
