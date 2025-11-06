#pragma once

#include <ModelTriangle.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include "Material.h"
#include "Vertex.h"
#include "Light.h"
#include <fstream>
#include <vector>
#include <unordered_map>
#include <Utils.h>
#include "RayTriangleIntersection.h"
#include "ExtraUtils.h"
#include "Ray.h"

#define MIN_DIST 0.001f
#define MAX_DIST 20.0f
#define AMBIENT 0.1f

class Scene {
public:
    Scene(const std::string &objFilename, std::vector<Light> &lights, float modelScale = 0.35);
    RayTriangleIntersection closestIntersection(Ray &ray, float minDist = MIN_DIST, float maxDist = MAX_DIST);
    glm::vec3 backgroundColour(float x, float y);
    void moveLight(Direction dir);
    std::vector<ModelTriangle> triangles{};
    std::unordered_map<std::string, Material> materials{};
    std::unordered_map<std::string, TextureMap> textures{};
    std::unordered_map<std::string, TextureMap> normalMaps{};
    std::vector<Light> &lights;
    TextureMap* background{nullptr};
private:
    void readObj(const std::string &objFilename, float modelScale);
    void readMtl(const std::string &filename);
};
