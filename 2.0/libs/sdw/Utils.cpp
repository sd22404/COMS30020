#include "Utils.h"

std::vector<std::string> split(const std::string &line, const char delimiter) {
	auto haystack = line;
	std::vector<std::string> tokens;
	size_t pos;
	while ((pos = haystack.find(delimiter)) != std::string::npos) {
		tokens.push_back(haystack.substr(0, pos));
		haystack.erase(0, pos + 1);
	}
	// Push the remaining chars onto the vector
	tokens.push_back(haystack);
	return tokens;
}

// Uses Cramer’s rule to convert from 2D coordinates to Barycentric proportional proximities
glm::vec3 convertToBarycentricCoordinates(const glm::vec2 v0, const glm::vec2 v1, const glm::vec2 v2, const glm::vec2 r)
{
    const glm::vec2 e0 = v1 - v0;
    const glm::vec2 e1 = v2 - v0;
    const glm::vec2 e2 = r - v0;
    const float d00 = glm::dot(e0, e0);
    const float d01 = glm::dot(e0, e1);
    const float d11 = glm::dot(e1, e1);
    const float d20 = glm::dot(e2, e0);
    const float d21 = glm::dot(e2, e1);
    const float denominator = d00 * d11 - d01 * d01;
    const float u = (d11 * d20 - d01 * d21) / denominator;
    const float v = (d00 * d21 - d01 * d20) / denominator;
    const float w = 1.0f - u - v;
    return {u,v,w};
}

float edgeFunction(const CanvasPoint &v0, const CanvasPoint &v1, const CanvasPoint &v2) {
    return (v2.x - v0.x) * (v1.y - v0.y) - (v2.y - v0.y) * (v1.x - v0.x);
}