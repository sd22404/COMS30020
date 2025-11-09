#pragma once

#include <ModelTriangle.h>
#include <TextureMap.h>
#include "Material.h"
#include "Light.h"
#include <unordered_map>
#include "HitInfo.h"
#include "ExtraUtils.h"
#include "Ray.h"
#include "Model.h"
#include "Obj.h"

#define MIN_DIST 0.001f
#define MAX_DIST 20.0f
#define AMBIENT 0.1f

class Scene {
public:
    Scene(const std::vector<Obj> &objs, std::vector<Light> &lights);
    HitInfo closestIntersection(const Ray &ray, float minDist = MIN_DIST, float maxDist = MAX_DIST) const;
    glm::vec3 backgroundColour(const int x, const int y) const;
    void moveLight(Direction dir) const;
    std::vector<Model> models{};
    std::vector<Light> &lights;
    const TextureMap *background{nullptr};
private:
    void readObj(const Obj &obj);
    std::unordered_map<std::string, Material> readMtl(const std::string &filename);
};
