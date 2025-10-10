#pragma once

#include <ModelTriangle.h>
#include <TextureMap.h>
#include "Material.h"
#include "Vertex.h"
#include <fstream>
#include <vector>
#include <unordered_map>
#include <Utils.h>

class Scene {
public:
    std::vector<ModelTriangle> triangles;
    std::unordered_map<std::string, Material> materials;
    std::unordered_map<std::string, TextureMap> textures;
    std::unordered_map<std::string, TextureMap> normalMaps;

    void readObj(const std::string &objFilename, float modelScale);
    void readMtl(const std::string &filename);
    Scene(const std::string &objFilename, float modelScale = 0.35);
};