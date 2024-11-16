#pragma once

#include <iostream>

struct TexturePoint {
	float x{};
	float y{};

	TexturePoint();
	TexturePoint(float xPos, float yPos);
	friend std::ostream &operator<<(std::ostream &os, const TexturePoint &point);
	friend TexturePoint operator*(const float &scale, const TexturePoint &point);
	friend TexturePoint operator+(const TexturePoint &other, const TexturePoint &point);
};
