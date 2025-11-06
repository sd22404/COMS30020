#pragma once

#include <DrawingWindow.h>
#include <Utils.h>
#include "Scene.h"
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <CanvasLine.h>
#include <Utils.h>
#include <algorithm>
#include "Camera.h"
#include "Ray.h"

class Renderer {
private:
    DrawingWindow &window;
    RenderMode rMode;
    LightingMode lMode;
    bool drawLight;
    std::vector<std::vector<float>> depthBuffer;
public:
    Renderer(DrawingWindow &window) :
        window(window), rMode(WIREFRAME), lMode(FLAT), drawLight(false), depthBuffer(std::vector<std::vector<float>>(window.height, std::vector<float>(window.width, 0.0f))) {}
    void setRenderMode(RenderMode mode) { rMode = mode; }
    void setLightingMode(LightingMode mode) { lMode = mode; }
    void toggleLight() { drawLight = !drawLight; }
    void draw(Scene &scene, Camera &cam);
private:
    void drawLine(CanvasPoint &p0, CanvasPoint &p1, uint32_t &colour);
    void drawTriangle(CanvasTriangle &triangle, uint32_t &colour);
    void fillTriangle(CanvasTriangle &triangle, uint32_t &colour);
    uint32_t traceRay(Ray &ray, Scene &scene, Camera &cam);
    void wireframe(Scene &scene, Camera &cam);
    void raster(Scene &scene, Camera &cam);
    void raytrace(Scene &scene, Camera &cam);
};
