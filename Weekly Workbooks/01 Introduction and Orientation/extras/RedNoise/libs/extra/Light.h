#pragma once

#include <ostream>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

#define DEFAULT_INT 25.0f
#define DEFAULT_LIGHT_POS glm::vec3(0.0f, 0.9f, 0.0f)

#define V0 glm::vec3(-0.64901096, 2.7384973, -0.51796794) * 0.35f - glm::vec3(0, 0.05, 0)
#define V1 glm::vec3(-0.64901096, 2.739334, 0.532032) * 0.35f - glm::vec3(0, 0.05, 0)
#define V2 glm::vec3(0.650989, 2.7384973, -0.51796794) * 0.35f - glm::vec3(0, 0.05, 0)

enum Type {POINT, AREA};

struct Light {
    Type type;
    glm::vec3 defaultPosition{DEFAULT_LIGHT_POS};
    glm::vec3 position{V0};
    glm::vec3 e1{V1 - V0};
    glm::vec3 e2{V2 - V0};
    float intensity{DEFAULT_INT};

    explicit Light(Type type);
    Light(Type type, float intensity);
    Light(Type type, glm::vec3 position, float intensity);
    Light(Type type, glm::vec3 v0, glm::vec3 e1, glm::vec3 e2, float intensity);

    glm::vec3 sample();

    friend std::ostream &operator<<(std::ostream &os, const Light &light);
};