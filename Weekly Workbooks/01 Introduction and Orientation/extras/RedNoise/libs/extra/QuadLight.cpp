#include "QuadLight.h"

QuadLight::QuadLight() : Light(DEFAULT_INT, V0) {}

QuadLight::QuadLight(glm::vec3 v0, glm::vec3 e1, glm::vec3 e2, float intensity) :
    e1(e1),
    e2(e2),
    Light(intensity, v0) {}
