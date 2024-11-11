#ifndef LIGHT_H
#define LIGHT_H
#include <ostream>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>


struct Light {
    glm::vec3 defaultPosition;
    float intensity;
    float radius;
    glm::vec3 position;

    Light();
    Light(float intensity, float radius);
    Light(float intensity, float radius, glm::vec3 position);

    friend std::ostream &operator<<(std::ostream &os, const Light &light);
};



#endif //LIGHT_H
