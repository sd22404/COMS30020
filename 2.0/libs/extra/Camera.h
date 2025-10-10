#pragma once

#include "CanvasPoint.h"
#include <DrawingWindow.h>
#include <glm/glm.hpp>

enum Direction {
    UP,
    DOWN,
    LEFT,
    RIGHT,
    FORWARD,
    BACKWARD
};

class Camera {
public:
    glm::vec3 position = glm::vec3(0, 0, 2);
    glm::vec3 rotation = glm::vec3(0, 0, 0);
    float focalLength = 1.0f;
    float speed = 0.1f;
    int width;
    int height;

    Camera(int width, int height) : width(width), height(height) {};
    CanvasPoint projectVertex(const glm::vec3 &vertex, float canvasScale);
    void move(Direction dir);
};