#pragma once

#include "ModelTriangle.h"
#include "TextureMap.h"
#include <ExtraUtils.h>

struct Model {
    std::string name;
    std::vector<ModelTriangle> triangles{};
    Material material{};

    glm::vec3 position{};
    glm::vec3 rotation{};
    glm::vec3 scale{1.0f, 1.0f, 1.0f};

    ShadingMode sMode = FLAT;
};