#include "Renderer.h"

void Renderer::draw(Scene &scene) {
	window.clearPixels();
    switch(mode) {
        case WIREFRAME:
            wireframe(scene);
            break;
        case RASTERISED:
    		raster(scene);
            break;
        case RAYTRACED:
    		raytrace(scene);
            break;
        default:
            break;
    }
}

void Renderer::drawLine(CanvasPoint &p0, CanvasPoint &p1, uint32_t &colour) {
    // number of steps is one more than the maximum difference in x- or y-direction
	float numberOfValues = glm::max(abs(p1.x - p0.x), abs(p1.y - p0.y));

	// interp along x and y, store in vector
	CanvasLine line = CanvasLine(p0, p1, numberOfValues);
	if (line.isOffScreen(window.width, window.height)) return;

	// set pixels based on interpolated values (rounded)
	for (int i = 0; i < numberOfValues; i++) {
		CanvasPoint &point = line.points[i];
		// check point is on screen
		if (!point.isOffScreen(window.width, window.height)) {
			window.setPixelColour(floor(point.x), floor(point.y), colour);
		}
	}
}

void Renderer::drawTriangle(CanvasTriangle &triangle, uint32_t &colour) {
    // draw lines between each pair of vertices
    drawLine(triangle.v0(), triangle.v1(), colour);
    drawLine(triangle.v1(), triangle.v2(), colour);
    drawLine(triangle.v2(), triangle.v0(), colour);
}

void Renderer::fillTriangle(CanvasTriangle &triangle, uint32_t &colour) {
	if (triangle.isOffScreen(window.width, window.height)) return;

	float maxX = glm::max(glm::max(triangle.v0().x, triangle.v1().x), triangle.v2().x);
	float maxY = glm::max(glm::max(triangle.v0().y, triangle.v1().y), triangle.v2().y);
	float minX = glm::min(glm::min(triangle.v0().x, triangle.v1().x), triangle.v2().x);
	float minY = glm::min(glm::min(triangle.v0().y, triangle.v1().y), triangle.v2().y);
	// fill in triangle
	for (float x = minX; x < maxX; x++) {
		for (float y = minY; y < maxY; y++) {
			glm::vec3 bCoords = convertToBarycentricCoordinates(
				{triangle.v0().x, triangle.v0().y},
				{triangle.v1().x, triangle.v1().y},
				{triangle.v2().x, triangle.v2().y},
				{x, y});
			float depth = bCoords.x * triangle.v1().depth + bCoords.y * triangle.v2().depth + bCoords.z * triangle.v0().depth;
			auto point = CanvasPoint(x, y, depth);
			if (!point.isOffScreen(window.width, window.height) &&
				bCoords.x >= 0 && bCoords.y >= 0 && bCoords.z >= 0 &&
				depth > depthBuffer[floor(y)][floor(x)] &&
				depth > 0) {
					window.setPixelColour(floor(x), floor(y), colour);
					depthBuffer[floor(y)][floor(x)] = depth;
			}
		}
	}
}

void Renderer::wireframe(Scene &scene) {
    for (auto &triangle : scene.triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = scene.cam.projectVertex(triangle.v0().position, window.height / 2.0f);
		CanvasPoint v1 = scene.cam.projectVertex(triangle.v1().position, window.height / 2.0f);
		CanvasPoint v2 = scene.cam.projectVertex(triangle.v2().position, window.height / 2.0f);
		CanvasTriangle canvasTriangle = {v0, v1, v2};
    	uint32_t colour =
			(255 << 24) +
			(int(triangle.material.colour.r * 255) << 16) +
			(int(triangle.material.colour.g * 255) << 8) +
			(int(triangle.material.colour.b * 255));
		drawTriangle(canvasTriangle, colour);
	}
}

void Renderer::raster(Scene &scene) {
	for(auto& row : depthBuffer) std::fill(row.begin(), row.end(), 0);
	for (auto &triangle : scene.triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = scene.cam.projectVertex(triangle.v0().position, window.height / 2.0f);
		CanvasPoint v1 = scene.cam.projectVertex(triangle.v1().position, window.height / 2.0f);
		CanvasPoint v2 = scene.cam.projectVertex(triangle.v2().position, window.height / 2.0f);
		CanvasTriangle canvasTriangle = {v0, v1, v2};
		uint32_t colour =
			(255 << 24) +
			(int(triangle.material.colour.r * 255) << 16) +
			(int(triangle.material.colour.g * 255) << 8) +
			(int(triangle.material.colour.b * 255));
		fillTriangle(canvasTriangle, colour);
	}
}

void Renderer::raytrace(Scene &scene) {
	for (float y = 0; y < window.height; y++) {
		for (float x = 0; x < window.width; x++) {
			glm::vec3 dir = scene.cam.projectRay(x, y, window.height / 2.0f);
			auto intersection = scene.closestIntersection(scene.cam.position, dir);
			auto triangle = intersection.intersectedTriangle;
			bool miss = intersection.triangleIndex == -1;
			uint32_t colour = 0;
			if (miss) colour = scene.backgroundColour(x, y);
			else colour =
				(255 << 24) +
				(int(triangle.material.colour.r * 255) << 16) +
				(int(triangle.material.colour.g * 255) << 8) +
				(int(triangle.material.colour.b * 255));
			window.setPixelColour(floor(x), floor(y), colour);
		}
	}
}
