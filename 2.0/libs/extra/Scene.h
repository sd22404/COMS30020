#pragma once

#include <ModelTriangle.h>
#include <TextureMap.h>
#include "Material.h"
#include "Vertex.h"
#include <fstream>
#include <vector>
#include <unordered_map>
#include <Utils.h>
#include "RayTriangleIntersection.h"
#include "Camera.h"

#define MIN_DIST 0.001f
#define MAX_DIST 20.0f

class Scene {
public:
    Scene(const std::string &objFilename, Camera &camera, float modelScale = 0.35);
    RayTriangleIntersection closestIntersection(glm::vec3 start, glm::vec3 dir, float minDist = MIN_DIST, float maxDist = MAX_DIST);
    uint32_t backgroundColour(float x, float y);
    std::vector<ModelTriangle> triangles{};
    std::unordered_map<std::string, Material> materials{};
    std::unordered_map<std::string, TextureMap> textures{};
    std::unordered_map<std::string, TextureMap> normalMaps{};
    Camera &cam;
    TextureMap* background{nullptr};
private:
    void readObj(const std::string &objFilename, float modelScale);
    void readMtl(const std::string &filename);
    static glm::vec3 vertexNormal(Vertex &vertex, std::vector<ModelTriangle> &triangles);
};
