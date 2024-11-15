#ifndef VERTEX_H
#define VERTEX_H
#include <ostream>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

struct Vertex {
    int index{};
    glm::vec3 position{};
    glm::vec3 normal{};

    Vertex();
    Vertex(float x, float y, float z);
    Vertex(float x, float y, float z, glm::vec3 normal);
    Vertex(glm::vec3 position);
    Vertex(glm::vec3 position, glm::vec3 normal);
    friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex);
};

#endif //VERTEX_H