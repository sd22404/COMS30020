#pragma once
#include <glm/glm.hpp>
#include <vector>

#define AREA_LIGHT_SAMPLES 4

#define CEIL0 (glm::vec3(-0.64901096, 2.7384973, -0.51796794) * 0.35f)
#define CEIL1 (glm::vec3(-0.64901096, 2.739334, 0.532032) * 0.35f)
#define CEIL2 (glm::vec3(0.650989, 2.7384973, -0.51796794) * 0.35f)

enum LightType {
    POINT,
    AREA
};

struct Light {
    Light(const glm::vec3 &position, const float &intensity)
        : position(position), colour(glm::vec3(1, 1, 1)), intensity(intensity) {
            sample();
        }
    Light(const glm::vec3 &position, const glm::vec3 &colour, const float &intensity)
        : position(position), colour(colour), intensity(intensity) {
            sample();
        }
    Light(const glm::vec3 &position, const float &intensity, const glm::vec3 &uVec, const glm::vec3 &vVec, const int uSteps = AREA_LIGHT_SAMPLES, const int vSteps = AREA_LIGHT_SAMPLES)
        : position(position), colour({1, 1, 1}), intensity(intensity), uVec(uVec), vVec(vVec), uSteps(uSteps), vSteps(vSteps), type(AREA) {
            sample();
        }
    Light(const glm::vec3 &position, const glm::vec3 &colour, const float &intensity, const glm::vec3 &uVec, const glm::vec3 &vVec, const int uSteps = AREA_LIGHT_SAMPLES, const int vSteps = AREA_LIGHT_SAMPLES)
        : position(position), colour(colour), intensity(intensity), uVec(uVec), vVec(vVec), uSteps(uSteps), vSteps(vSteps), type(AREA) {
            sample();
        }

    glm::vec3 position;
    glm::vec3 colour;
    float intensity;

    LightType type{POINT};
    glm::vec3 uVec{0, 0, 0};
    glm::vec3 vVec{0, 0, 0};
    int uSteps{1};
    int vSteps{1};

    std::vector<glm::vec3> samples;
    void sample();
};