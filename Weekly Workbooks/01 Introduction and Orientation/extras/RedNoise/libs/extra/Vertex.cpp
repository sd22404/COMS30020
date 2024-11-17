#include "Vertex.h"

Vertex::Vertex() = default;

Vertex::Vertex(float x, float y, float z) : index(-1), position(glm::vec3(x, y, z)), normal(), texturePoint() {}

Vertex::Vertex(float x, float y, float z, glm::vec3 normal) : index(-1), position(glm::vec3(x, y, z)), normal(normal), texturePoint() {}

Vertex::Vertex(float x, float y, float z, glm::vec3 normal, TexturePoint tp) : index(-1), position(glm::vec3(x, y, z)), normal(normal), texturePoint(tp) {}

Vertex::Vertex(glm::vec3 position) : index(-1), position(position), normal(), texturePoint() {}

Vertex::Vertex(glm::vec3 position, glm::vec3 normal) : index(-1), position(position), normal(normal), texturePoint() {}

Vertex::Vertex(glm::vec3 position, glm::vec3 normal, TexturePoint tp) : index(-1), position(position), normal(normal), texturePoint(tp) {}

std::ostream &operator<<(std::ostream &os, const Vertex &vertex) {
    os << to_string(vertex.position) << ", " << to_string(vertex.normal);
    return os;
}
