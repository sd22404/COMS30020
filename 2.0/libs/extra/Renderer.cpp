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
    drawLine(triangle[0], triangle[1], colour);
    drawLine(triangle[1], triangle[2], colour);
    drawLine(triangle[2], triangle[0], colour);
}

void Renderer::fillTriangle(CanvasTriangle &triangle, uint32_t &colour, const TextureMap &textureMap) {
	if (triangle.isOffScreen(window.width, window.height)) return;

	float area = edgeFunction(triangle[0], triangle[1], triangle[2]);

	int maxX = ceil(std::max({triangle[0].x, triangle[1].x, triangle[2].x}));
	int maxY = ceil(std::max({triangle[0].y, triangle[1].y, triangle[2].y}));
	int minX = floor(std::min({triangle[0].x, triangle[1].x, triangle[2].x}));
	int minY = floor(std::min({triangle[0].y, triangle[1].y, triangle[2].y}));
	// fill in triangle
	for (int x = minX; x < maxX; x++) {
		for (int y = minY; y < maxY; y++) {

			CanvasPoint p = CanvasPoint(x, y);

			const float u = edgeFunction(triangle[1], triangle[2], p) / area;
            const float v = edgeFunction(triangle[2], triangle[0], p) / area;
            const float w = 1.0f - u - v;

			float depth = u * triangle[0].depth + v * triangle[1].depth + w * triangle[2].depth;
			p.depth = depth;
			if (!p.isOffScreen(window.width, window.height) &&
				u >= 0 && v >= 0 && w >= 0 &&
				depth > depthBuffer[y][x]) {
					if (!textureMap.name.empty()) {
						int texX = (int)(floor(u * triangle[0].texturePoint.x + v * triangle[1].texturePoint.x + w * triangle[2].texturePoint.x));
						int texY = (int)(floor(u * triangle[0].texturePoint.y + v * triangle[1].texturePoint.y + w * triangle[2].texturePoint.y));
						window.setPixelColour(x, y, textureMap.pixels[texY * textureMap.width + texX]);
					} else {
						window.setPixelColour(x, y, colour);
					}
					depthBuffer[y][x] = depth;
			}
		}
	}
}

glm::vec3 Renderer::shade(RayTriangleIntersection &hit, Light &light, glm::vec3 &N, glm::vec3 &V) {
	ModelTriangle &triangle = hit.intersectedTriangle;
	const Material &mat = triangle.material;
	const float u = hit.u, v = hit.v;

	glm::vec3 lightDir = light.position - hit.intersectionPoint;
	float dist = glm::length(lightDir);

	glm::vec3 L = glm::normalize(lightDir);

	float attenuation = light.intensity / (4 * M_PIf * dist * dist + 1.0f);
	float diffuse;
	float specular;

	if (lMode == GOURAUD) {
		float diff0 = std::max(dot(triangle[0].normal, L), 0.0f);
		float diff1 = std::max(dot(triangle[1].normal, L), 0.0f);
		float diff2 = std::max(dot(triangle[2].normal, L), 0.0f);
		diffuse = (u * diff1 + v * diff2 + (1.0f - u - v) * diff0);

		float spec0 = std::pow(std::max(dot(glm::reflect(-L, triangle[0].normal), V), 0.0f), mat.shininess);
		float spec1 = std::pow(std::max(dot(glm::reflect(-L, triangle[1].normal), V), 0.0f), mat.shininess);
		float spec2 = std::pow(std::max(dot(glm::reflect(-L, triangle[2].normal), V), 0.0f), mat.shininess);
		specular = (u * spec1 + v * spec2 + (1.0f - u - v) * spec0);

		return diffuse * mat.diffuse + specular * mat.specular * light.colour * attenuation;
	}

	diffuse = std::max(dot(N, L), 0.0f);

	glm::vec3 R = glm::reflect(-L, N);
	specular = std::pow(std::max(dot(R, V), 0.0f), mat.shininess);

	return diffuse * mat.diffuse + specular * mat.specular * light.colour * attenuation;
}

glm::vec3 Renderer::traceRay(Ray &ray, Scene &scene, int depth) {
	auto hit = scene.closestIntersection(ray);
	bool miss = hit.triangleIndex == -1;
	if (miss) return scene.backgroundColour(0, 0);
	ModelTriangle &triangle = hit.intersectedTriangle;
	const Material &mat = triangle.material;
	glm::vec3 surface = hit.intersectionPoint;
	glm::vec3 colour = mat.diffuse * mat.ambient;
	glm::vec3 V = -normalize(ray.dir);
	glm::vec3 N = (lMode == PHONG) ?
		normalize(hit.u * triangle[1].normal + hit.v * triangle[2].normal + (1.0f - hit.u - hit.v) * triangle[0].normal)
		: triangle.normal;

	for (auto light: scene.lights) {
		glm::vec3 lightDir = light.position - surface;
		float dist = glm::length(lightDir);

		Ray shadowRay(surface, normalize(lightDir), dist);
		auto shadowHit = scene.closestIntersection(shadowRay, MIN_DIST, dist - MIN_DIST);
		bool inShadow = shadowHit.triangleIndex != -1;
		if (inShadow) continue;

		colour += shade(hit, light, N, V);
	}

	if (depth < MAX_DEPTH && mat.reflectivity > 0.0f) {
		glm::vec3 R = glm::reflect(-V, N);
		Ray reflectRay(surface + 0.001f * R, R);
		glm::vec3 reflectColor = traceRay(reflectRay, scene, depth + 1);
		colour = colour * (1.0f - mat.reflectivity) + reflectColor * mat.reflectivity;
	}

	colour = glm::clamp(colour, glm::vec3(0.0f), glm::vec3(1.0f));
	return colour;
}


void Renderer::wireframe(Scene &scene, Camera &cam) {
    for (auto &triangle : scene.triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = cam.projectVertex(triangle[0].position, window.height / 2.0f);
		CanvasPoint v1 = cam.projectVertex(triangle[1].position, window.height / 2.0f);
		CanvasPoint v2 = cam.projectVertex(triangle[2].position, window.height / 2.0f);
		CanvasTriangle canvasTriangle = {v0, v1, v2};
    	uint32_t colour = packColour(triangle.material.diffuse);
		drawTriangle(canvasTriangle, colour);
	}
}

void Renderer::raster(Scene &scene, Camera &cam) {
	for(auto& row : depthBuffer) std::fill(row.begin(), row.end(), 0.0f);
	for (auto &triangle : scene.triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = cam.projectVertex(triangle[0], window.height / 2.0f);
		CanvasPoint v1 = cam.projectVertex(triangle[1], window.height / 2.0f);
		CanvasPoint v2 = cam.projectVertex(triangle[2], window.height / 2.0f);
		CanvasTriangle canvasTriangle = {v0, v1, v2};
		uint32_t colour = packColour(triangle.material.diffuse);
		if (!triangle.texture.empty()) {
			fillTriangle(canvasTriangle, colour, scene.textures[triangle.texture]);
		}
		else fillTriangle(canvasTriangle, colour);
	}
}

void Renderer::raytrace(Scene &scene, Camera &cam) {
	for (int y = 0; y < window.height; y++) {
		for (int x = 0; x < window.width; x++) {
			Ray ray = cam.projectRay(x, y, window.height / 2.0f);
			uint32_t colour = packColour(traceRay(ray, scene));
			window.setPixelColour(x, y, colour);
		}
	}
}
