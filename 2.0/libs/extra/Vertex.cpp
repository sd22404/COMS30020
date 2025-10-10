#include "Vertex.h"

Vertex::Vertex() = default;

Vertex::Vertex(float x, float y, float z) : index(-1), position(glm::vec3(x, y, z)) {}

Vertex::Vertex(float x, float y, float z, glm::vec3 normal) : index(-1), position(glm::vec3(x, y, z)), normal(normal) {}

Vertex::Vertex(float x, float y, float z, glm::vec3 normal, TexturePoint tp) : index(-1), position(glm::vec3(x, y, z)), normal(normal), texturePoint(tp) {}

Vertex::Vertex(glm::vec3 position) : index(-1), position(position) {}

Vertex::Vertex(glm::vec3 position, glm::vec3 normal) : index(-1), position(position), normal(normal) {}

Vertex::Vertex(glm::vec3 position, glm::vec3 normal, TexturePoint tp) : index(-1), position(position), normal(normal), texturePoint(tp) {}

Vertex::Vertex(glm::vec3 position, glm::vec3 normal, TexturePoint tp, TexturePoint np) : index(-1), position(position), normal(normal), texturePoint(tp), normalPoint(np) {}

std::ostream &operator<<(std::ostream &os, const Vertex &vertex) {
    os << glm::to_string(vertex.position) << ", " << glm::to_string(vertex.normal);
    return os;
}