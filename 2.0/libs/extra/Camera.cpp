#include "Camera.h"

CanvasPoint Camera::projectVertex(const Vertex &vertex, float canvasScale) {
    // vertex in terms of camera coordinates
	glm::vec3 finalPos = (vertex.position - position) * rotation;
	// transform onto image plane
	float u = -canvasScale * focalLength * finalPos.x / finalPos.z + width / 2.0f;
	float v = canvasScale * focalLength * finalPos.y / finalPos.z + height / 2.0f;
    CanvasPoint point(u, v, -1 / finalPos.z);
    point.texturePoint = vertex.texturePoint;
	return point;
}

Ray Camera::projectRay(int &x, int &y, float canvasScale) {
    // convert from SDL coordinate system into 3D/model coordinate system
    float sdlx = x - width / 2.0f;
    float sdlY = -y + height / 2.0f;
    // generate canvas point in 3D space, adjusted by cameraOrientation
    glm::vec3 canvasPoint3D = position + glm::vec3(sdlx, sdlY, -(focalLength * canvasScale)) * inverse(rotation);
    // subtract cameraPosition and normalise to get ray direction
    Ray ray(position, normalize(canvasPoint3D - position));
    return ray;
}

glm::mat3 Camera::rotateY(float angle) {
    return {
        glm::vec3(cos(angle), 0, -sin(angle)),
        glm::vec3(0, 1, 0),
        glm::vec3(sin(angle), 0, cos(angle))
    };
}

float Camera::degToRad(float deg) {
    return M_PI * deg / 180;
}

void Camera::lookAt(glm::vec3 target) {
    glm::vec3 forward = normalize(position - target);
    glm::vec3 right = normalize(cross(glm::vec3(0, 1, 0), forward));
    glm::vec3 up = normalize(cross(forward, right));
    rotation = glm::mat3(right, up, forward);
}

void Camera::move(Direction dir) {
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
