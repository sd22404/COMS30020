#include "Light.h"

Light::Light() : position(defaultPosition) {}

Light::Light(float intensity) :
    position(defaultPosition),
    intensity(intensity) {}

Light::Light(float intensity, glm::vec3 position) :
    position(position),
    intensity(intensity) {}

std::ostream &operator<<(std::ostream &os, const Light &light) {
    os << "(" << light.intensity << ") " << to_string(light.position);
    return os;
}
