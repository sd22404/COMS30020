#pragma once

#include "CanvasPoint.h"
#include <iostream>
#include <array>

struct CanvasTriangle {
	std::array<CanvasPoint, 3> vertices{};

	CanvasTriangle();
	CanvasTriangle(const CanvasPoint &v0, const CanvasPoint &v1, const CanvasPoint &v2);
	CanvasPoint operator[](size_t i) const;
	CanvasPoint &operator[](size_t i);
	bool isOffScreen(float width, float height);
	friend std::ostream &operator<<(std::ostream &os, const CanvasTriangle &triangle);
};

std::ostream &operator<<(std::ostream &os, const CanvasTriangle &triangle);
