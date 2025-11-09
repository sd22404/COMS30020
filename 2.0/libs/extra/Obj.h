#pragma once

#include "ExtraUtils.h"
#include <string>

struct Obj {
    std::string filename;
    ShadingMode shadingMode = FLAT;
    glm::vec3 offset = glm::vec3(0);
    float scale = 0.35f;
};