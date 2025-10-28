#include "TexturePoint.h"

TexturePoint::TexturePoint() = default;
TexturePoint::TexturePoint(float xPos, float yPos) : x(xPos), y(yPos) {}

std::ostream &operator<<(std::ostream &os, const TexturePoint &point) {
	os << "x: " << point.x << " y: " << point.y;
	return os;
}

TexturePoint operator*(const float &scale, const TexturePoint &point) {
	TexturePoint tp = TexturePoint(point.x * scale, point.y * scale);
	return tp;
}

TexturePoint operator+(const TexturePoint &other, const TexturePoint &point) {
	TexturePoint tp = TexturePoint(other.x + point.x, other.y + point.y);
	return tp;
}

TexturePoint operator-(const TexturePoint &other, const TexturePoint &point) {
	TexturePoint tp = TexturePoint(other.x - point.x, other.y - point.y);
	return tp;
}

bool operator!=(const TexturePoint &other, const TexturePoint &point) {
	return !(other.x == point.x && other.y == point.y);
}
