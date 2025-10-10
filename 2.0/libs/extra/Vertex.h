#pragma once

#include <ostream>
#include "../sdw/TexturePoint.h"
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

struct Vertex {
    int index{};
    glm::vec3 position{};
    glm::vec3 normal{};
    TexturePoint texturePoint{};
    TexturePoint normalPoint{};

    Vertex();
    Vertex(float x, float y, float z);
    Vertex(float x, float y, float z, glm::vec3 normal);
    Vertex(float x, float y, float z, glm::vec3 normal, TexturePoint texturePoint);
    explicit Vertex(glm::vec3 position);
    Vertex(glm::vec3 position, glm::vec3 normal);
    Vertex(glm::vec3 position, glm::vec3 normal, TexturePoint texturePoint);
    Vertex(glm::vec3 position, glm::vec3 normal, TexturePoint texturePoint, TexturePoint normalPoint);
    friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex);
};