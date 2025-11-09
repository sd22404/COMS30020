#pragma once

#include <glm/glm.hpp>
#include <iostream>

struct HitInfo {
	glm::vec3 intersectionPoint;
	float distanceFromCamera{};
	float u{}, v{}, w{};
	int modelIndex{-1};
	int triIndex{-1};

	HitInfo();
	HitInfo(const glm::vec3 &point, float distance, int modelIndex, int triangleIndex);
	friend std::ostream &operator<<(std::ostream &os, const HitInfo &intersection);
};
