#include "CanvasPoint.h"

bool CanvasPoint::isOffScreen(const float width, const float height) const {
	if (x >= width || x < 0 || y >= height || y < 0 || depth < 0) return true;
	return false;
}

CanvasPoint::CanvasPoint() :
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(const float xPos, const float yPos) :
		x(xPos),
		y(yPos),
		depth(0.0),
		brightness(1.0),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(const float xPos, const float yPos, const float pointDepth) :
		x(xPos),
		y(yPos),
		depth(pointDepth),
		brightness(1.0),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(const float xPos, const float yPos, const float pointDepth, const float pointBrightness) :
		x(xPos),
		y(yPos),
		depth(pointDepth),
		brightness(pointBrightness),
		texturePoint(-1, -1) {}

std::ostream &operator<<(std::ostream &os, const CanvasPoint &point) {
	os << "(" << point.x << ", " << point.y << ", " << point.depth << ") " << point.brightness;
	return os;
}
