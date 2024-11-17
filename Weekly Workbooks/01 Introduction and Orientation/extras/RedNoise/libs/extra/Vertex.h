#pragma once

#include <ostream>
#include <TexturePoint.h>
#include <glm/detail/type_vec3.hpp>
#include <glm/gtx/string_cast.hpp>

struct Vertex {
    int index{};
    glm::vec3 position{};
    glm::vec3 normal{};
    TexturePoint texturePoint{};

    Vertex();
    Vertex(float x, float y, float z);
    Vertex(float x, float y, float z, glm::vec3 normal);
    Vertex(float x, float y, float z, glm::vec3 normal, TexturePoint texturePoint);
    Vertex(glm::vec3 position);
    Vertex(glm::vec3 position, glm::vec3 normal);
    Vertex(glm::vec3 position, glm::vec3 normal, TexturePoint texturePoint);
    friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex);
};