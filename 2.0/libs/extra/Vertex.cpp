#include "Vertex.h"

Vertex::Vertex() = default;

Vertex::Vertex(const float x, const float y, const float z) : index(-1), position(glm::vec3(x, y, z)) {}

Vertex::Vertex(const float x, const float y, const float z, const glm::vec3 normal) : index(-1), position(glm::vec3(x, y, z)), normal(normal) {}

Vertex::Vertex(const float x, const float y, const float z, const glm::vec3 normal, const TexturePoint tp) : index(-1), position(glm::vec3(x, y, z)), normal(normal), texturePoint(tp) {}

Vertex::Vertex(const glm::vec3 position) : index(-1), position(position) {}

Vertex::Vertex(const glm::vec3 position, const glm::vec3 normal) : index(-1), position(position), normal(normal) {}

Vertex::Vertex(const glm::vec3 position, const glm::vec3 normal, const TexturePoint tp) : index(-1), position(position), normal(normal), texturePoint(tp) {}

std::ostream &operator<<(std::ostream &os, const Vertex &vertex) {
    os << glm::to_string(vertex.position) << ", " << glm::to_string(vertex.normal);
    return os;
}