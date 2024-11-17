#include "PointLight.h"

#define DEFAULT_INT 1.0f
#define DEFAULT_POS glm::vec3(0.0f, 0.9f, 0.0f)

PointLight::PointLight() :
    defaultPosition(DEFAULT_POS),
    position(defaultPosition),
    intensity(DEFAULT_INT) {}

PointLight::PointLight(float intensity) :
    defaultPosition(DEFAULT_POS),
    position(defaultPosition),
    intensity(intensity) {}

PointLight::PointLight(float intensity, glm::vec3 position) :
    defaultPosition(DEFAULT_POS),
    position(position),
    intensity(intensity) {}

std::ostream &operator<<(std::ostream &os, const PointLight &light) {
    os << "(" << light.intensity << ") " << to_string(light.position);
    return os;
}
