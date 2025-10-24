#pragma once

#include <DrawingWindow.h>
#include <Utils.h>
#include "Scene.h"
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <CanvasLine.h>
#include <Utils.h>
#include <algorithm>

class Renderer {
public:
    enum RenderMode {
        WIREFRAME,
        RASTERISED,
        RAYTRACED
    };
private:
    DrawingWindow &window;
    RenderMode mode;
    std::vector<std::vector<float>> depthBuffer;
public:
    Renderer(DrawingWindow &window) :
        window(window), mode(WIREFRAME), depthBuffer(std::vector<std::vector<float>>(window.height, std::vector<float>(window.width, 0.0f))) {}
    void setMode(RenderMode rMode) { mode = rMode; }
    void draw(Scene &scene);
private:
    void drawLine(CanvasPoint &p0, CanvasPoint &p1, uint32_t &colour);
    void drawTriangle(CanvasTriangle &triangle, uint32_t &colour);
    void fillTriangle(CanvasTriangle &triangle, uint32_t &colour);
    void wireframe(Scene &scene);
    void raster(Scene &scene);
    void raytrace(Scene &scene);
};
