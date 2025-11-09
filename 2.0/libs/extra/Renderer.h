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

#define MAX_DEPTH 4

class Renderer {
private:
    DrawingWindow &window;
    RenderMode rMode;
    bool drawLight;
    std::vector<std::vector<float>> depthBuffer;
public:
    explicit Renderer(DrawingWindow &window) :
        window(window), rMode(WIREFRAME), drawLight(false), depthBuffer(std::vector<std::vector<float>>(window.height, std::vector<float>(window.width, 0.0f))) {}
    void setRenderMode(const RenderMode mode) { rMode = mode; }
    void toggleLight() { drawLight = !drawLight; }
    void draw(const Scene &scene, const Camera &cam);
private:
    void drawLine(const CanvasPoint &p0, const CanvasPoint &p1, const uint32_t &colour) const;
    void drawTriangle(const CanvasTriangle &triangle, const uint32_t &colour) const;
    void fillTriangle(const CanvasTriangle &triangle, const uint32_t &colour, const TextureMap &textureMap = TextureMap());
    glm::vec3 shade(const HitInfo &hit, const ModelTriangle &triangle, const Model &model, const glm::vec3 &baseColour, const Light &light, const glm::vec3 &N, const glm::vec3 &V) const;
    glm::vec3 traceRay(const Ray &ray, const Scene &scene, int depth = 0);
    void wireframe(const Scene &scene, const Camera &cam) const;
    void raster(const Scene &scene, const Camera &cam);
    void raytrace(const Scene &scene, const Camera &cam);
};
