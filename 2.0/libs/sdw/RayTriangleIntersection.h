#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include "ModelTriangle.h"

struct RayTriangleIntersection {
	glm::vec3 intersectionPoint;
	float distanceFromCamera{};
	ModelTriangle intersectedTriangle;
	size_t triangleIndex{};
	float u{};
	float v{};
	uint32_t textureSample{};
	uint32_t normalSample{};

	RayTriangleIntersection();
	RayTriangleIntersection(const glm::vec3 &point, float distance, ModelTriangle triangle, size_t index);
	friend std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection);
};
