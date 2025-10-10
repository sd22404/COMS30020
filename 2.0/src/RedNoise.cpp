#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>

#define WIDTH 640
#define HEIGHT 480


void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			auto baryCoords = convertToBarycentricCoordinates(glm::vec2(0, HEIGHT), glm::vec2(WIDTH, HEIGHT), glm::vec2(WIDTH / 2, 0), glm::vec2(x, y));
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

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
