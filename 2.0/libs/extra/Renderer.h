#pragma once

#include <DrawingWindow.h>
#include <Utils.h>
#include "Scene.h"
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include "RenderUtils.h"
#include "Camera.h"

enum RenderMode {
    WIREFRAME,
    RASTERISED,
    RAYTRACED
};

class Renderer {
private:
    DrawingWindow &window;
    RenderMode mode;
    Scene &scene;
    Camera &camera;
public:
    Renderer(DrawingWindow &window, Scene &scene, Camera &cam) : window(window), mode(WIREFRAME), scene(scene), camera(cam) {}
    void draw(Scene &scene);
    void drawLine(CanvasPoint p0, CanvasPoint p1, glm::vec3 colour);
    void drawTriangle(CanvasTriangle triangle, glm::vec3 colour);
    void wireframe();
};
