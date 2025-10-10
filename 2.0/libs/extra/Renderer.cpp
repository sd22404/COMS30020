#include "Renderer.h"

void Renderer::draw(Scene &scene) {
    switch(mode) {
        case WIREFRAME:
            wireframe();
            break;
        case RASTERISED:
            break;
        case RAYTRACED:
            break;
        default:
            break;
    }
}

void Renderer::drawLine(CanvasPoint p0, CanvasPoint p1, glm::vec3 colour) {
    // CULL COMPLETELY OFF-SCREEN LINES
	// if (isOffScreen(p0, p1)) return;

    // number of steps is one more than the maximum difference in x- or y-direction
	int numberOfValues = ceil(glm::max(abs(p1.x - p0.x), abs(p1.y - p0.y)) + 1);

	// interp along x and y, store in vector
	std::vector<CanvasPoint> line = RenderUtils::interpolateCanvasPoints(p0, p1, numberOfValues);
	// pack colour
	uint32_t packedCol = (255 << 24) + (int(colour.r * 255) << 16) + (int(colour.g * 255) << 8) + int(colour.b * 255);

	// set pixels based on interpolated values (rounded)
	for (int i = 0; i < numberOfValues; i++) {
		CanvasPoint point = CanvasPoint(line[i].x, line[i].y, line[i].depth);
		// check point is on screen
		// if (!isOffScreen(point)) {
			window.setPixelColour(floor(point.x), floor(point.y), packedCol);
		// }
	}
}

void Renderer::drawTriangle(CanvasTriangle triangle, glm::vec3 colour) {
    // draw lines between each pair of vertices
    drawLine(triangle.v0(), triangle.v1(), colour);
    drawLine(triangle.v1(), triangle.v2(), colour);
    drawLine(triangle.v2(), triangle.v0(), colour);
}

void Renderer::wireframe() {
    window.clearPixels();
    for (auto triangle: scene.triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = camera.projectVertex(triangle.v0().position, window.height / 2.0f);
		CanvasPoint v1 = camera.projectVertex(triangle.v1().position, window.height / 2.0f);
		CanvasPoint v2 = camera.projectVertex(triangle.v2().position, window.height / 2.0f);
		CanvasTriangle canvasTriangle = {v0, v1, v2};
		drawTriangle(canvasTriangle, triangle.material.colour);
	}
}

/*
void Renderer::colourfulTriangle() {
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
*/
