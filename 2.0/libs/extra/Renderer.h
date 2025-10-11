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
    std::vector<std::vector<float>> depthBuffer;
public:
    Renderer(DrawingWindow &window, Scene &scene, Camera &cam) :
        window(window), mode(WIREFRAME), scene(scene), camera(cam), depthBuffer(std::vector<std::vector<float>>(window.height, std::vector<float>(window.width, 0))) {}
    void setMode(RenderMode rMode) { mode = rMode; }
    void draw(Scene &scene);
private:
    void drawLine(CanvasPoint &p0, CanvasPoint &p1, uint32_t &colour);
    void drawTriangle(CanvasTriangle &triangle, uint32_t &colour);
    void fillTriangle(CanvasTriangle &triangle, uint32_t &colour);
    void wireframe();
    void raster();
};
