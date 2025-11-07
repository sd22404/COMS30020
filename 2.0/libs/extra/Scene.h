#pragma once

#include <ModelTriangle.h>
#include <TextureMap.h>
#include "Material.h"
#include "Light.h"
#include <vector>
#include <unordered_map>
#include "RayTriangleIntersection.h"
#include "ExtraUtils.h"
#include "Ray.h"

#define MIN_DIST 0.001f
#define MAX_DIST 20.0f
#define AMBIENT 0.1f

class Scene {
public:
    Scene(const std::string &objFilename, std::vector<Light> &lights, float modelScale = 0.35);
    RayTriangleIntersection closestIntersection(const Ray &ray, float minDist = MIN_DIST, float maxDist = MAX_DIST) const;
    glm::vec3 backgroundColour(float x, float y) const;
    void moveLight(Direction dir) const;
    std::vector<ModelTriangle> triangles{};
    std::unordered_map<std::string, Material> materials{};
    std::unordered_map<std::string, TextureMap> textures{};
    std::unordered_map<std::string, TextureMap> normalMaps{};
    std::vector<Light> &lights;
    const TextureMap *background{nullptr};
private:
    void readObj(const std::string &objFilename, float modelScale);
    void readMtl(const std::string &filename);
};
