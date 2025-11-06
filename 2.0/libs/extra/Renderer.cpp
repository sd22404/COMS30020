#include "Renderer.h"

void Renderer::draw(Scene &scene, Camera &cam) {
	window.clearPixels();
    switch(rMode) {
        case WIREFRAME:
            wireframe(scene, cam);
            break;
        case RASTERISED:
    		raster(scene, cam);
            break;
        case RAYTRACED:
    		raytrace(scene, cam);
            break;
        default:
            break;
    }

	if (drawLight) {
		CanvasPoint ceilingLight = cam.projectVertex(scene.lights.at(0).position, window.height / 2.0f);
		CanvasPoint left = CanvasPoint(ceilingLight.x - 5, ceilingLight.y - 5);
		CanvasPoint right = CanvasPoint(ceilingLight.x + 5, ceilingLight.y - 5);
		CanvasTriangle lightTriangle = {ceilingLight, left, right};
		uint32_t colour = 0xFFFF0000;
		drawTriangle(lightTriangle, colour);
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

	float area = edgeFunction(triangle.v0(), triangle.v1(), triangle.v2());

	int maxX = ceil(std::max({triangle.v0().x, triangle.v1().x, triangle.v2().x}));
	int maxY = ceil(std::max({triangle.v0().y, triangle.v1().y, triangle.v2().y}));
	int minX = floor(std::min({triangle.v0().x, triangle.v1().x, triangle.v2().x}));
	int minY = floor(std::min({triangle.v0().y, triangle.v1().y, triangle.v2().y}));
	// fill in triangle
	for (int x = minX; x < maxX; x++) {
		for (int y = minY; y < maxY; y++) {

			CanvasPoint p = CanvasPoint(x, y);

			const float u = edgeFunction(triangle.v1(), triangle.v2(), p) / area;
            const float v = edgeFunction(triangle.v2(), triangle.v0(), p) / area;
            const float w = 1.0f - u - v;

			float depth = u * triangle.v0().depth + v * triangle.v1().depth + w * triangle.v2().depth;
			p.depth = depth;
			if (!p.isOffScreen(window.width, window.height) &&
				u >= 0 && v >= 0 && w >= 0 &&
				depth > depthBuffer[y][x]) {
					window.setPixelColour(x, y, colour);
					depthBuffer[y][x] = depth;
			}
		}
	}
}

static float proximityBrightness(const float distToLight, const float intensity) {
	return intensity / (4 * M_PIf * distToLight * distToLight + 5);
}

static float angularBrightness(const glm::vec3 &normal, const glm::vec3 &toLight) {
	return dot(normal, normalize(toLight));
}

static float specularBrightness(const glm::vec3 &normal, const glm::vec3 &toLight, const glm::vec3 &toCam, const float n = 32) {
	glm::vec3 view = normalize(toCam);
	glm::vec3 incidence = normalize(toLight);
	glm::vec3 reflection = incidence - 2.0f * normal * dot(incidence, normal);
	return glm::pow(glm::max(dot(reflection, view), 0.0f), n);
}

static float flatBrightness(const glm::vec3 &normal, float intensity, const glm::vec3 &toLight, const glm::vec3 &toCam) {
	float brightness = proximityBrightness(glm::length(toLight), intensity);
	brightness *= angularBrightness(normal, toLight);
	brightness += specularBrightness(normal, toLight, toCam);
	return brightness;
}

static float gouraudBrightness(ModelTriangle &triangle, const float u, const float v, const float intensity, const glm::vec3 &toLight, const glm::vec3 &toCam) {
	glm::vec3 &v0Norm = triangle.v0().normal;
	glm::vec3 &v1Norm = triangle.v1().normal;
	glm::vec3 &v2Norm = triangle.v2().normal;
	float v0Brightness = flatBrightness(v0Norm, intensity, toLight, toCam);
	float v1Brightness = flatBrightness(v1Norm, intensity, toLight, toCam);
	float v2Brightness = flatBrightness(v2Norm, intensity, toLight, toCam);
	float brightness = u * v1Brightness + v * v2Brightness + (1.0f - u - v) * v0Brightness;
	return brightness;
}

static float phongBrightness(ModelTriangle &triangle, const float u, const float v, const float intensity, const glm::vec3 &toLight, const glm::vec3 &toCam) {
	glm::vec3 normal = normalize(u * triangle.v1().normal + v * triangle.v2().normal + (1.0f - u - v) * triangle.v0().normal);
	return flatBrightness(normal, intensity, toLight, toCam);
}

uint32_t Renderer::traceRay(Ray &ray, Scene &scene, Camera &cam) {
	auto intersection = scene.closestIntersection(ray);
	bool miss = intersection.triangleIndex == -1;
	if (miss) return scene.backgroundColour(0, 0);
	ModelTriangle &triangle = intersection.intersectedTriangle;
	glm::vec3 surface = intersection.intersectionPoint;
	glm::vec3 toCam = cam.position - surface;
	float totalBrightness = 0;

	for (auto light: scene.lights) {
		float brightness = 0;
		glm::vec3 toLight = light.position - surface;

		switch (lMode) {
			case FLAT:
				brightness = flatBrightness(triangle.normal, light.intensity, toLight, toCam);
				break;
			case GOURAUD: {
				brightness = gouraudBrightness(triangle, intersection.u, intersection.v, light.intensity, toLight, toCam);
				break;
			}
			case PHONG:
				brightness = phongBrightness(triangle, intersection.u, intersection.v, light.intensity, toLight, toCam);
				break;
			default:
				break;
		}

		//shadow
		Ray shadowRay(surface, normalize(toLight), glm::length(toLight));
		auto shadowHit = scene.closestIntersection(shadowRay, MIN_DIST, shadowRay.dist - MIN_DIST);
		bool inShadow = shadowHit.triangleIndex != -1;
		if (inShadow) continue;

		totalBrightness += brightness;
	}

	totalBrightness = glm::clamp(totalBrightness, AMBIENT, 1.0f);

	return (255 << 24) + (
		   (int(totalBrightness * triangle.material.colour.r * 255) << 16) +
		   (int(totalBrightness * triangle.material.colour.g * 255) << 8) +
		   (int(totalBrightness * triangle.material.colour.b * 255)));
}


void Renderer::wireframe(Scene &scene, Camera &cam) {
    for (auto &triangle : scene.triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = cam.projectVertex(triangle.v0().position, window.height / 2.0f);
		CanvasPoint v1 = cam.projectVertex(triangle.v1().position, window.height / 2.0f);
		CanvasPoint v2 = cam.projectVertex(triangle.v2().position, window.height / 2.0f);
		CanvasTriangle canvasTriangle = {v0, v1, v2};
    	uint32_t colour =
			(255 << 24) +
			(int(triangle.material.colour.r * 255) << 16) +
			(int(triangle.material.colour.g * 255) << 8) +
			(int(triangle.material.colour.b * 255));
		drawTriangle(canvasTriangle, colour);
	}
}

void Renderer::raster(Scene &scene, Camera &cam) {
	for(auto& row : depthBuffer) std::fill(row.begin(), row.end(), 0.0f);
	for (auto &triangle : scene.triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = cam.projectVertex(triangle.v0().position, window.height / 2.0f);
		CanvasPoint v1 = cam.projectVertex(triangle.v1().position, window.height / 2.0f);
		CanvasPoint v2 = cam.projectVertex(triangle.v2().position, window.height / 2.0f);
		CanvasTriangle canvasTriangle = {v0, v1, v2};
		uint32_t colour =
			(255 << 24) +
			(int(triangle.material.colour.r * 255) << 16) +
			(int(triangle.material.colour.g * 255) << 8) +
			(int(triangle.material.colour.b * 255));
		fillTriangle(canvasTriangle, colour);
	}
}

void Renderer::raytrace(Scene &scene, Camera &cam) {
	for (int y = 0; y < window.height; y++) {
		for (int x = 0; x < window.width; x++) {
			Ray ray = cam.projectRay(x, y, window.height / 2.0f);
			uint32_t colour = traceRay(ray, scene, cam);
			window.setPixelColour(x, y, colour);
		}
	}
}
