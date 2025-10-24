#pragma once

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "CanvasPoint.h"

std::vector<std::string> split(const std::string &line, char delimiter);
glm::vec3 convertToBarycentricCoordinates(glm::vec2 v0, glm::vec2 v1, glm::vec2 v2, glm::vec2 r);
float edgeFunction(const CanvasPoint &v0, const CanvasPoint &v1, const CanvasPoint &v2);