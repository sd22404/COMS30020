#pragma once

#include "CanvasPoint.h"
#include <glm/glm.hpp>
#include "glm/gtx/string_cast.hpp"
#include "Enums.h"
#include "Ray.h"

class Camera {
private:
    int width;
    int height;
    bool orbiting = false;
    glm::vec3 startPosition = glm::vec3(0, 0, 4);
public:
    glm::vec3 position = glm::vec3(0, 0, 4);
    glm::mat3 rotation = glm::mat3(glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec3(0, 0, 1));
    float focalLength = 2.0f;
    float speed = 0.2f;

    Camera(int width, int height, float focalLength, glm::vec3 position) : width(width), height(height), startPosition(position), position(position), focalLength(focalLength) {};
    CanvasPoint projectVertex(const glm::vec3 &vertex, float canvasScale);
    Ray projectRay(int &x, int &y, float canvasScale);
    void reset();
    void lookAt(glm::vec3 target);
    void move(Direction dir);
    void toggleOrbit();
    void orbit();

    static glm::mat3 rotateY(float angle);
    static float degToRad(float deg);
};