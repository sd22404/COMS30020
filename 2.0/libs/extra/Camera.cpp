#include "Camera.h"

CanvasPoint Camera::projectVertex(const glm::vec3 &vertex, float canvasScale) {
    // vertex in terms of camera coordinates
	glm::vec3 vertexFromCamera = vertex - position;
	// adjust by cameraOrientation
	glm::vec3 adjustedPos = vertexFromCamera * rotation;
	float x = adjustedPos.x; float y = adjustedPos.y; float z = -adjustedPos.z; // invert z so that depth is positive
	// transform onto image plane
	float u = canvasScale * focalLength * x / abs(z) + width / 2.0f; // abs z to keep points that are behind camera
	float v = -canvasScale * focalLength * y / abs(z) + height / 2.0f;
	return {u, v, 1/z};
};

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
};