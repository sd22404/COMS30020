#pragma once

#include <ostream>
#include <glm/detail/type_vec3.hpp>

#include "Light.h"

#define V0 glm::vec3(-0.64901096, 2.7384973, -0.51796794) * 0.35f
#define V1 glm::vec3(-0.64901096, 2.739334, 0.532032) * 0.35f
#define V2 glm::vec3(0.650989, 2.7384973, -0.51796794) * 0.35f

struct QuadLight : Light {
    glm::vec3 e1{V1 - V0};
    glm::vec3 e2{V2 - V0};

    QuadLight();
    QuadLight(glm::vec3 v0, glm::vec3 e1, glm::vec3 e2, float intensity);
};
