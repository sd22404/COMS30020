#pragma once

#include <ostream>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

#define DEFAULT_INT 10.0f
#define DEFAULT_LIGHT_POS glm::vec3(0.0f, 0.9f, 0.0f)

#define CEIL0 (glm::vec3(-0.64901096, 2.7384973, -0.51796794) * 0.35f)
#define CEIL1 (glm::vec3(-0.64901096, 2.739334, 0.532032) * 0.35f)
#define CEIL2 (glm::vec3(0.650989, 2.7384973, -0.51796794) * 0.35f)

#define IDENT0 (glm::vec3(2.864524, -0.507219, 0.650000) * 0.35f - glm::vec3(0.05, 0, 0))
#define IDENT1 (glm::vec3(3.135476, 0.507219, 0.650000) * 0.35f - glm::vec3(0.05, 0, 0))
#define IDENT2 (glm::vec3(2.864524, -0.507219, -0.650000) * 0.35f - glm::vec3(0.05, 0, 0))

enum Type {POINT, AREA};

struct Light {
    Type type;
    glm::vec3 defaultPosition{DEFAULT_LIGHT_POS};
    glm::vec3 position{CEIL0};
    glm::vec3 e1{CEIL1 - CEIL0};
    glm::vec3 e2{CEIL2 - CEIL0};
    float intensity{DEFAULT_INT};

    explicit Light(Type type);
    Light(Type type, float intensity);
    Light(Type type, glm::vec3 position, float intensity);
    Light(Type type, glm::vec3 v0, glm::vec3 e1, glm::vec3 e2, float intensity);

    glm::vec3 sample();

    friend std::ostream &operator<<(std::ostream &os, const Light &light);
};