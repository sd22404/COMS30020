#include "Light.h"

#define DEFAULT_INT 1.0f
#define DEFAULT_RAD 1.0f
#define DEFAULT_POS glm::vec3(0.0f, 0.9f, 0.0f)
#define DEFAULT_AMBIENT 0.2f

Light::Light() :
    defaultPosition(DEFAULT_POS),
    intensity(DEFAULT_INT),
    radius(DEFAULT_RAD),
    ambient(DEFAULT_AMBIENT),
    position(defaultPosition) {}

Light::Light(float intensity, float radius) :
    defaultPosition(DEFAULT_POS),
    intensity(intensity),
    radius(radius),
    ambient(DEFAULT_AMBIENT),
    position(defaultPosition) {}

Light::Light(float intensity, float radius, glm::vec3 position) :
    defaultPosition(DEFAULT_POS),
    intensity(intensity),
    radius(radius),
    ambient(DEFAULT_AMBIENT),
    position(position) {}

std::ostream &operator<<(std::ostream &os, const Light &light) {
    os << "(" << light.intensity << ", " << light.radius << ") " << to_string(light.position);
    return os;
}
