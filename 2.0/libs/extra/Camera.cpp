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
    const glm::vec3 dir = glm::vec3(sdlX, sdlY, -(focalLength * canvasScale)) * inverse(rotation);
    // subtract cameraPosition and normalise to get ray direction
    return {position, normalize(dir)};
}

static glm::mat3 rotateAround(const glm::vec3 &axis, const float theta) {
    glm::vec3 u = normalize(axis);
	glm::mat3 W = {glm::vec3(0.0f, -u.z, u.y), glm::vec3(u.z, 0.0f, -u.x), glm::vec3(-u.y, u.x, 0.0f)};
	glm::mat3 I = {glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)};
	return I + sinf(theta) * W + 2.0f * powf(sinf(theta / 2.0f), 2.0f) * W * W;
}

static float degToRad(const float deg) {
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
            position += speed * normalize(rotation[1]);
            break;
        case DOWN:
            position -= speed * normalize(rotation[1]);
            break;
        case LEFT:
            position -= speed * normalize(rotation[0]);
            break;
        case RIGHT:
            position += speed * normalize(rotation[0]);
            break;
        case FORWARD:
            position -= speed * normalize(rotation[2]);
            break;
        case BACKWARD:
            position += speed * normalize(rotation[2]);
            break;
        default:
            break;
    }
}

void Camera::rotate(const Direction dir) {
    float theta = degToRad(speed * 10.0f);
    switch (dir) {
        case LEFT:
            rotation = rotateAround({0, 1, 0}, -theta) * rotation;
            break;
        case RIGHT:
            rotation = rotateAround({0, 1, 0}, theta) * rotation;
            break;
        case UP:
            rotation = rotateAround(rotation[0], -theta) * rotation;
            break;
        case DOWN:
            rotation = rotateAround(rotation[0], theta) * rotation;
            break;
        default:
            return;
    }
}

void Camera::reset() {
    position = startPosition;
    rotation = startRotation;
}

void Camera::toggleOrbit() { orbiting = !orbiting; }

void Camera::orbit() {
    if (!orbiting) return;
    position = rotateAround({0, 1, 0}, degToRad(speed * 10.0f)) * position;
    lookAt({0, 0, 0});
}
