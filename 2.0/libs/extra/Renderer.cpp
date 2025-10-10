#include "Renderer.h"

void Renderer::draw() {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			auto baryCoords = convertToBarycentricCoordinates(glm::vec2(0, window.height), glm::vec2(window.width, window.height), glm::vec2(window.width / 2, 0), glm::vec2(x, y));
			float red, green, blue;
			if (baryCoords.x < 0 || baryCoords.y < 0 || baryCoords.z < 0) {
				red = 0; green = 0; blue = 0;
			} else {
				red = glm::clamp(baryCoords.x * 255.0f, 0.0f, 255.0f);
				green = glm::clamp(baryCoords.y * 255.0f, 0.0f, 255.0f);
				blue = glm::clamp(baryCoords.z * 255.0f, 0.0f, 255.0f);
			}
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}