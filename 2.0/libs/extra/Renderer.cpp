#include "Renderer.h"

void Renderer::draw(const Scene &scene, const Camera &cam) {
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
		const auto ceilingLight = cam.projectVertex(Vertex(scene.lights.at(0).position), static_cast<float>(window.height) / 2.0f);
		const auto left = CanvasPoint(ceilingLight.x - 5, ceilingLight.y - 5);
		const auto right = CanvasPoint(ceilingLight.x + 5, ceilingLight.y - 5);
		const CanvasTriangle lightTriangle = {ceilingLight, left, right};
		constexpr uint32_t colour = 0xFFFF0000;
		drawTriangle(lightTriangle, colour);
	}
}

void Renderer::drawLine(const CanvasPoint &p0, const CanvasPoint &p1, const uint32_t &colour) const {
    // number of steps is one more than the maximum difference in x- or y-direction
	const long numberOfValues = std::lroundf(glm::max(abs(p1.x - p0.x), abs(p1.y - p0.y)));

	// interpolate along x and y, store in vector
	auto line = CanvasLine(p0, p1, numberOfValues);
	if (line.isOffScreen(static_cast<float>(window.width), static_cast<float>(window.height))) return;

	// set pixels based on interpolated values (rounded)
	for (int i = 0; i < numberOfValues; i++) {
		CanvasPoint &point = line.points[i];
		// check point is on screen
		if (!point.isOffScreen(static_cast<float>(window.width), static_cast<float>(window.height))) {
			window.setPixelColour(floor(point.x), floor(point.y), colour);
		}
	}
}

void Renderer::drawTriangle(const CanvasTriangle &triangle, const uint32_t &colour) const {
    // draw lines between each pair of vertices
    drawLine(triangle[0], triangle[1], colour);
    drawLine(triangle[1], triangle[2], colour);
    drawLine(triangle[2], triangle[0], colour);
}

void Renderer::fillTriangle(const CanvasTriangle &triangle, const uint32_t &colour, const TextureMap &textureMap) {
	// if (triangle.isOffScreen(static_cast<float>(window.width), static_cast<float>(window.height))) return;

	const float area = edgeFunction(triangle[0], triangle[1], triangle[2]);

	const int maxX = ceil(std::max({triangle[0].x, triangle[1].x, triangle[2].x}));
	const int maxY = ceil(std::max({triangle[0].y, triangle[1].y, triangle[2].y}));
	const int minX = floor(std::min({triangle[0].x, triangle[1].x, triangle[2].x}));
	const int minY = floor(std::min({triangle[0].y, triangle[1].y, triangle[2].y}));
	// fill in triangle
	for (int x = minX; x < maxX; x++) {
		for (int y = minY; y < maxY; y++) {

			auto p = CanvasPoint(static_cast<float>(x), static_cast<float>(y));

			const float u = edgeFunction(triangle[1], triangle[2], p) / area;
            const float v = edgeFunction(triangle[2], triangle[0], p) / area;
            const float w = 1.0f - u - v;

			const float depth = u * triangle[0].depth + v * triangle[1].depth + w * triangle[2].depth;
			p.depth = depth;
			if (!p.isOffScreen(static_cast<float>(window.width), static_cast<float>(window.height)) &&
				u >= 0 && v >= 0 && w >= 0 &&
				depth > depthBuffer[y][x]) {
				if (!textureMap.name.empty()) {
					const float texX = u * triangle[0].texturePoint.x + v * triangle[1].texturePoint.x + w * triangle[2].texturePoint.x;
					const float texY = u * triangle[0].texturePoint.y + v * triangle[1].texturePoint.y + w * triangle[2].texturePoint.y;
					const int iTexX = std::floor(texX * static_cast<float>(textureMap.width - 1));
					const int iTexY = std::floor(texY * static_cast<float>(textureMap.height - 1));
					window.setPixelColour(x, y, textureMap.pixels[iTexY * textureMap.width + iTexX]);
				} else {
					window.setPixelColour(x, y, colour);
				}
				depthBuffer[y][x] = depth;
			}
		}
	}
}

glm::vec3 Renderer::shade(const HitInfo &hit, const ModelTriangle &triangle, const Model &model, const glm::vec3 &baseColour, const Light &light, const glm::vec3 &N, const glm::vec3 &V) const {
	const glm::vec3 lightDir = light.position - hit.intersectionPoint;
	const float dist = glm::length(lightDir);

	const glm::vec3 L = glm::normalize(lightDir);

	const float attenuation = light.intensity / (4 * M_PIf * dist * dist + 1.0f);
	float diffuse;
	float specular;

	if (model.sMode == GOURAUD) {
		const float diff0 = std::max(dot(triangle[0].normal, L), 0.0f);
		const float diff1 = std::max(dot(triangle[1].normal, L), 0.0f);
		const float diff2 = std::max(dot(triangle[2].normal, L), 0.0f);
		diffuse = (hit.u * diff1 + hit.v * diff2 + hit.w * diff0);

		const float spec0 = std::pow(std::max(dot(glm::reflect(-L, triangle[0].normal), V), 0.0f), model.material.shininess);
		const float spec1 = std::pow(std::max(dot(glm::reflect(-L, triangle[1].normal), V), 0.0f), model.material.shininess);
		const float spec2 = std::pow(std::max(dot(glm::reflect(-L, triangle[2].normal), V), 0.0f), model.material.shininess);
		specular = (hit.u * spec1 + hit.v * spec2 + hit.w * spec0);

		return (diffuse * baseColour + specular * model.material.specular) * light.colour * attenuation;
	}

	diffuse = std::max(dot(N, L), 0.0f);

	const glm::vec3 R = glm::reflect(-L, N);
	specular = std::pow(std::max(dot(R, V), 0.0f), model.material.shininess);

	return (diffuse * baseColour + specular * model.material.specular) * light.colour * attenuation;
}

glm::vec3 Renderer::traceRay(const Ray &ray, const Scene &scene, int depth) {
	auto hit = scene.closestIntersection(ray);
	if (hit.modelIndex == -1) return scene.backgroundColour(0, 0);

	const Model &model = scene.models[hit.modelIndex];
	const ModelTriangle &triangle = model.triangles[hit.triIndex];
	const Material &mat = model.material;
	glm::vec3 P = hit.intersectionPoint;

	const float texX = hit.u * triangle[1].texturePoint.x + hit.v * triangle[2].texturePoint.x + hit.w * triangle[0].texturePoint.x;
	const float texY = hit.u * triangle[1].texturePoint.y + hit.v * triangle[2].texturePoint.y + hit.w * triangle[0].texturePoint.y;

	// base colour (diffuse or texture sample)
	glm::vec3 baseColour = mat.diffuse;
	if (!mat.texture.name.empty()) {
		const size_t iTexX = std::floor(texX * static_cast<float>(mat.texture.width));
		const size_t iTexY = std::floor(texY * static_cast<float>(mat.texture.height));
		const size_t idx = (iTexY * mat.texture.width + iTexX) % (mat.texture.pixels.size());
		const uint32_t textureSample = mat.texture.pixels[idx];
		baseColour = unpackColour(textureSample);
	}

	// compute normal at intersection
	glm::vec3 V = -normalize(ray.dir);
	glm::vec3 N = (model.sMode == PHONG)
		? normalize(hit.u * triangle[1].normal + hit.v * triangle[2].normal + hit.w * triangle[0].normal)
		: triangle.normal;
	if (!mat.normalMap.name.empty()) {
		const size_t iTexX = std::floor(texX * static_cast<float>(mat.normalMap.width));
		const size_t iTexY = std::floor(texY * static_cast<float>(mat.normalMap.height));
		const size_t idx = (iTexY * mat.normalMap.width + iTexX) % (mat.normalMap.pixels.size());
		const uint32_t normalSample = mat.normalMap.pixels[idx];
		glm::vec3 N_ts = unpackColour(normalSample) * 2.0f - glm::vec3(1.0f);
		N = normalize(N + N_ts);
	}

	glm::vec3 colour = baseColour * mat.ambient;
	for (const auto &light : scene.lights) {
		glm::vec3 Ldir = light.position - P;
		float dist = glm::length(Ldir);
		
		if (model.sMode == FLAT) {
			Ray shadowRay(P + N * 0.001f, normalize(Ldir), dist);
			auto shadowHit = scene.closestIntersection(shadowRay, MIN_DIST, dist - MIN_DIST);
			if (shadowHit.modelIndex != -1) continue;
		}
		
		colour += shade(hit, triangle, model, baseColour, light, N, V);
	}

	if (depth < MAX_DEPTH && mat.reflectivity > 0.0f) {
		glm::vec3 R = glm::reflect(-V, N);
		Ray reflectRay(P + 0.001f * R, R);
		glm::vec3 reflectColor = traceRay(reflectRay, scene, depth + 1);
		colour = colour * (1.0f - mat.reflectivity) + reflectColor * mat.reflectivity;
	}

	return glm::clamp(colour, glm::vec3(0.0f), glm::vec3(1.0f));
}


void Renderer::wireframe(const Scene &scene, const Camera &cam) const {
	for (auto &model : scene.models) {
		for (auto &triangle : model.triangles) {
			// for each model triangle, project vertices onto canvas and draw resulting triangle
			CanvasPoint v0 = cam.projectVertex(triangle[0], static_cast<float>(window.height) / 2.0f);
			CanvasPoint v1 = cam.projectVertex(triangle[1], static_cast<float>(window.height) / 2.0f);
			CanvasPoint v2 = cam.projectVertex(triangle[2], static_cast<float>(window.height) / 2.0f);
			CanvasTriangle canvasTriangle = {v0, v1, v2};
			uint32_t colour = packColour(model.material.diffuse);
			drawTriangle(canvasTriangle, colour);
		}
	}
}

void Renderer::raster(const Scene &scene, const Camera &cam) {
	for (auto& row : depthBuffer) std::fill(row.begin(), row.end(), 0.0f);
	for (auto &model : scene.models) {
		for (auto &triangle : model.triangles) {
			// for each model triangle, project vertices onto canvas and draw resulting triangle
			CanvasPoint v0 = cam.projectVertex(triangle[0], static_cast<float>(window.height) / 2.0f);
			CanvasPoint v1 = cam.projectVertex(triangle[1], static_cast<float>(window.height) / 2.0f);
			CanvasPoint v2 = cam.projectVertex(triangle[2], static_cast<float>(window.height) / 2.0f);
			CanvasTriangle canvasTriangle = {v0, v1, v2};
			uint32_t colour = packColour(model.material.diffuse);
			!model.material.texture.name.empty()
			? fillTriangle(canvasTriangle, colour, model.material.texture)
			: fillTriangle(canvasTriangle, colour);
		}
	}
}

void Renderer::raytrace(const Scene &scene, const Camera &cam) {
	for (int y = 0; y < static_cast<int>(window.height); y++) {
		for (int x = 0; x < static_cast<int>(window.width); x++) {
			Ray ray = cam.projectRay(x, y, static_cast<float>(window.height) / 2.0f);
			const uint32_t colour = packColour(traceRay(ray, scene));
			window.setPixelColour(x, y, colour);
		}
	}
}
