#include "Light.h"

Light::Light(Type type) : type(type) {}

Light::Light(Type type, float intensity) : type(type), intensity(intensity) {}

Light::Light(Type type, glm::vec3 position, float intensity) : type(type), intensity(intensity), position(position) {}

Light::Light(Type type, glm::vec3 v0, glm::vec3 e1, glm::vec3 e2, float intensity) : type(type), intensity(intensity), position(v0), e1(e1), e2(e2) {}

glm::vec3 Light::sample() {
    glm::vec3 sample;
    switch (type) {
        case POINT:
            {
                sample = position;
                break;
            }
        case AREA:
            {
                float u = float(rand() % 101) / 100.0f, v = float(rand() % 101) / 100.0f;
                sample = position + u * e1 + v * e2;
                break;
            }
        default:
            break;
    }

    return sample;
}


std::ostream &operator<<(std::ostream &os, const Light &light) {
    os << "(type: " << light.type << ") " << to_string(light.position) << ", " << light.intensity;
    return os;
}
