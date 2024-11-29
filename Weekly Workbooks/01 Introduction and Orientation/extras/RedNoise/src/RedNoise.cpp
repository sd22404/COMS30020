#include <algorithm>
#include <map>
#include <vector>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Material.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <TextureMap.h>
#include <glm/glm.hpp>
#include <Camera.h>
#include <functional>
#include <Light.h>
#include <sstream>

#include "omp.h"

#define WIDTH 640
#define HEIGHT 480

#define MIN_INTERSECT_DISTANCE 0.001f
#define MAX_INTERSECT_DISTANCE 20.0f
#define AREA_LIGHT_SAMPLES 2
#define MAX_RAY_BOUNCES 4

const std::string SCENE_ONE = "./assets/scenes/scene-1.obj";
const std::string SCENE_TWO = "./assets/scenes/scene-2.obj";

std::function<float(RayTriangleIntersection &surface, glm::vec3 &lightPos, float intensity, Camera &cam)> brightnessFunction;

RayTriangleIntersection closestIntersection(std::vector<ModelTriangle> &triangles, glm::vec3 &startPoint, glm::vec3 &rayDirection, Camera &cam, int bounces, float maxDist = MAX_INTERSECT_DISTANCE, float minDist = MIN_INTERSECT_DISTANCE) {
	float inverseClosestDistance = 0;
	// initialise empty RayTriangleIntersection
	RayTriangleIntersection intersection;
	intersection.triangleIndex = -1;
	// for each triangle, check for possible solution to intersection equation
	for (size_t i = 0; i < triangles.size(); i++) {
		ModelTriangle &triangle = triangles[i];
		// calculate edge vectors
		glm::vec3 e0 = triangle.v1().position - triangle.v0().position;
		glm::vec3 e1 = triangle.v2().position - triangle.v0().position;
		// calculate vector from startPoint to triangle
		glm::vec3 SPVector = startPoint - triangle.v0().position;
		// generate direction/edge matrix
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		// find possible solution in [t, u, v]
		glm::vec3 possibleSolution = inverse(DEMatrix) * SPVector;
		float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
		// if closer than previously found solution, and within the bounds of the triangle, set new closest intersection
		if (t > minDist && t < maxDist && 1 / t > inverseClosestDistance && u >= 0 && u <= 1.0 && v >= 0 && v <= 1.0 && (u + v) <= 1.0) {
			intersection = RayTriangleIntersection(startPoint + t * rayDirection, t, triangle, i, {u, v, 1.0f - u - v});
			inverseClosestDistance = 1 / t;
		}
	}

	// interpolate normal if phong shading
	if (cam.shadingMode == PHONG) {
		ModelTriangle &triangle = intersection.intersectedTriangle;
		glm::vec3 &bCoords = intersection.barycentricPoint;
		intersection.intersectedTriangle.normal = normalize(bCoords.x * triangle.v1().normal + bCoords.y * triangle.v2().normal + bCoords.z * triangle.v0().normal);
	}

	// bounce again if mirrored surface
	ModelTriangle &triangle = intersection.intersectedTriangle;
	if (triangle.material.mirrored && bounces > 0) {
		glm::vec3 normal = triangle.normal, incidence = rayDirection, start = intersection.intersectionPoint;
		glm::vec3 reflect = normalize(incidence - 2.0f * normal * dot(incidence, normal));
		return closestIntersection(triangles, start, reflect, cam, bounces - 1);
	}
	if (triangle.material.glassy && bounces > 0) {
		glm::vec3 normal = triangle.normal, incidence = rayDirection, start = intersection.intersectionPoint, newDirection;
		glm::vec3 reflect = normalize(incidence - 2.0f * normal * dot(incidence, normal));

		// determine which side of triangle and adjust ri
		float ri = dot(incidence, normal) > 0 ? 1.0f / triangle.material.refractiveIndex : triangle.material.refractiveIndex;
		// using Snell's law
		float cosTheta = glm::min(dot(-incidence, normal), 1.0f);
		float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
		glm::vec3 refractPerpendicular = ri * (incidence + cosTheta * normal);
		glm::vec3 refractParallel = -sqrt(float(abs(1.0f - pow(refractPerpendicular.length(), 2)))) * normal;
		glm::vec3 refract = normalize(refractPerpendicular + refractParallel);
		// using Schlick approximation for Fresnel
		float r0 = powf((1 - ri) / (1 + ri), 2);
		float reflectance = r0 + (1 - r0) * powf(1 - cosTheta, 5);
		// check for total internal reflection
		if (ri * sinTheta > 1.0f || reflectance > 1.0f) newDirection = reflect; // reflectance > rand() / (float) RAND_MAX
		else newDirection = refract;

		return closestIntersection(triangles, start, newDirection, cam, bounces - 1);
	}

	return intersection;
}

glm::vec3 rayFromCanvasPoint(CanvasPoint &canvasPoint, Camera &cam, float canvasScale = HEIGHT / 2.0) {
	// convert from SDL coordinate system into 3D/model coordinate system
	float x = canvasPoint.x - WIDTH / 2.0f;
	float y = -canvasPoint.y + HEIGHT / 2.0f;
	// generate canvas point in 3D space, adjusted by cameraOrientation
	glm::vec3 canvasPoint3D = cam.position + glm::vec3(x, y, -(cam.focalLength * canvasScale)) * inverse(cam.rotation);
	// subtract cameraPosition and normalise to get ray direction
	glm::vec3 ray = normalize(canvasPoint3D - cam.position);

	return ray;
}

bool inShadow(std::vector<ModelTriangle> &triangles, glm::vec3 &surface, glm::vec3 &lightPos, Camera &cam) {
	// calculate shadowRay direction
	glm::vec3 shadowRay = lightPos - surface;
	// get the closest intersection of shadowRay from surface
	RayTriangleIntersection obstacle = closestIntersection(triangles, surface, shadowRay, cam, MAX_RAY_BOUNCES, length(shadowRay));
	// if obstacle doesn't exist then no shadow and isn't a light itself
	if (obstacle.triangleIndex == size_t(-1) || obstacle.intersectedTriangle.material.emissive) return false;
	return true;
}

float proximityLighting(glm::vec3 &point, glm::vec3 &lightPos, float intensity) {
	float dist = powf(distance(point, lightPos), 2);
	// + constant in denominator to avoid infinite brightness at light source
	float brightness = intensity / float(4 * M_PI * dist + 5);
	return brightness;
}

float angularLighting(glm::vec3 &point, glm::vec3 &normal, glm::vec3 &lightPos) {
	float brightness = dot(normal, normalize(lightPos - point));
	return brightness;
}

float specularLighting(glm::vec3 &point, glm::vec3 &normal, glm::vec3 &lightPos, glm::vec3 &camPos, float n = 256) {
	glm::vec3 view = normalize(camPos - point);
	glm::vec3 incidence = normalize(point - lightPos);
	glm::vec3 reflection = incidence - 2.0f * normal * dot(incidence, normal);
	float brightness = pow(glm::max(dot(reflection, view), 0.0f), n);
	return brightness;
}

float surfaceBrightness(glm::vec3 &point, glm::vec3 &normal, glm::vec3 &lightPos, float intensity, Camera &cam) {
	// calculate proximity lighting
	float prox = proximityLighting(point, lightPos, intensity);
	// calculate angle of incidence lighting
	float ang = angularLighting(point, normal, lightPos);
	// calculate specular lighting
	float spec = specularLighting(point, normal, lightPos, cam.position);

	float brightness = prox * ang + spec;
	brightness = glm::clamp(brightness, 0.0f, 1.0f);
	return brightness;
}

float gouraudBrightness(RayTriangleIntersection &intersection, glm::vec3 &lightPos, float intensity, Camera &cam) {
	glm::vec3 &bCoords = intersection.barycentricPoint;
	glm::vec3 &v0_norm = intersection.intersectedTriangle.v0().normal;
	glm::vec3 &v1_norm = intersection.intersectedTriangle.v1().normal;
	glm::vec3 &v2_norm = intersection.intersectedTriangle.v2().normal;
	float v0_brightness = surfaceBrightness(intersection.intersectionPoint, v0_norm, lightPos, intensity, cam);
	float v1_brightness = surfaceBrightness(intersection.intersectionPoint, v1_norm, lightPos, intensity, cam);
	float v2_brightness = surfaceBrightness(intersection.intersectionPoint, v2_norm, lightPos, intensity, cam);
	float brightness = bCoords.x * v1_brightness + bCoords.y * v2_brightness + bCoords.z * v0_brightness;
	return brightness;
}

float flatBrightness(RayTriangleIntersection &intersection, glm::vec3 &lightPos, float intensity, Camera &cam) {
	return surfaceBrightness(intersection.intersectionPoint, intersection.intersectedTriangle.normal, lightPos, intensity, cam);
}


TexturePoint barycentricTexturePoint(ModelTriangle &triangle, glm::vec3 &bCoords) {
	TexturePoint tp = bCoords.x * triangle.v1().texturePoint + bCoords.y * triangle.v2().texturePoint + bCoords.z * triangle.v0().texturePoint;
	return tp;
}

TexturePoint barycentricNormalPoint(ModelTriangle &triangle, glm::vec3 &bCoords) {
	TexturePoint np = bCoords.x * triangle.v1().normalPoint + bCoords.y * triangle.v2().normalPoint + bCoords.z * triangle.v0().normalPoint;
	return np;
}


CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 &vertexPosition, Camera &cam, float canvasScale = HEIGHT / 2.0) {
	// vertex in terms of camera coordinates
	glm::vec3 vertexFromCamera = vertexPosition - cam.position;
	// adjust by cameraOrientation
	glm::vec3 adjustedPos = vertexFromCamera * cam.rotation;
	float x = adjustedPos.x; float y = adjustedPos.y; float z = -adjustedPos.z; // invert z so that depth is positive
	// transform onto image plane
	float u = canvasScale * cam.focalLength * x / abs(z) + WIDTH / 2.0f; // abs z to keep points that are behind camera
	float v = -canvasScale * cam.focalLength * y / abs(z) + HEIGHT / 2.0f;
	return {u, v, 1/z};
}

// return rotation matrix in each axis of input theta
glm::mat3 rotateX(float theta) {
	return {glm::vec3(1, 0, 0), glm::vec3(0, cos(theta), sin(theta)), glm::vec3(0, -sin(theta), cos(theta))};
}
glm::mat3 rotateY(float theta) {
	return {glm::vec3(cos(theta), 0, -sin(theta)), glm::vec3(0, 1, 0), glm::vec3(sin(theta), 0, cos(theta))};
}
glm::mat3 rotateZ(float theta) {
	return {glm::vec3(cos(theta), sin(theta), 0), glm::vec3(-sin(theta), cos(theta), 0), glm::vec3(0, 0, 1)};
}
glm::mat3 rotateV(glm::vec3 vector, float theta) {
	glm::vec3 u = normalize(vector);
	glm::mat3 W = {glm::vec3(0, -u.z, u.y), glm::vec3(u.z, 0, -u.x), glm::vec3(-u.y, u.x, 0)};
	glm::mat3 I = {glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec3(0, 0, 1)};
	return I + sin(theta) * W + 2.0f * powf(sin(theta / 2.0f), 2) * W * W;
}

float degToRad(float theta) {
	return theta * 2.0f * float(M_PI) / 180.0f;
}

void lookAt(glm::vec3 &point, Camera &cam) {
	// set camera axes to align with difference in current position and given point
	glm::vec3 forward = normalize(cam.position - point);
	glm::vec3 right = normalize(cross(glm::vec3(0, 1, 0), forward));
	glm::vec3 up = normalize(cross(forward, right));
	cam.rotation = glm::mat3(right, up, forward);
}

bool isOffScreen(CanvasPoint &v0, const CanvasPoint &v1 = CanvasPoint(), const CanvasPoint &v2 = CanvasPoint()) {
	// checks null case against brightness - will need changing if points with zero brightness are used
	if ((v0.x >= WIDTH && (v1.brightness == 0 || v1.x >= WIDTH) && (v2.brightness == 0 || v2.x >= WIDTH)) ||
		(v0.x < 0 && (v1.brightness == 0 || v1.x < 0) && (v2.brightness == 0 || v2.x < 0)) ||
		(v0.y >= HEIGHT && (v1.brightness == 0 || v1.y >= HEIGHT) && (v2.brightness == 0 || v2.y >= HEIGHT)) ||
		(v0.y < 0 && (v1.brightness == 0 || v1.y < 0) && (v2.brightness == 0 || v2.y < 0)) ||
		// if depth is negative, point is behind camera
		(v0.depth < 0 && (v1.brightness == 0 || v1.depth < 0) && (v2.brightness == 0 || v2.depth < 0))) return true;
	return false;
}


std::vector<float> interpolateSingleFloats(float from, float to, size_t numberOfValues) {
	// if one or fewer values, return only the start point (as it will be the same as the end point)
	std::vector<float> result;
	if (numberOfValues < 2) return result;

	float interval = to - from;
	float step = interval / float(numberOfValues - 1);

	for (size_t i = 1; i < numberOfValues; i++) {
		result.push_back(from + float(i) * step);
	}
	return result;
}

std::vector<TexturePoint> interpolateTexturePoints(TexturePoint from, TexturePoint to, size_t numberOfValues) {
	// if one or fewer values, return only the start point (as it will be the same as the end point)
	std::vector<TexturePoint> result = {from};
	if (numberOfValues < 2) return result;

	TexturePoint interval = {to.x - from.x, to.y - from.y};
	float stepX = interval.x / float(numberOfValues - 1);
	float stepY = interval.y / float(numberOfValues - 1);

	for (size_t i = 1; i < numberOfValues; i++) {
		TexturePoint point = {from.x + float(i) * stepX, from.y + float(i) * stepY};
		result.push_back(point);
	}
	return result;
}

std::vector<CanvasPoint> interpolateCanvasPoints(CanvasPoint from, CanvasPoint to, size_t numberOfValues) {
	// if one or fewer values, return only the start point (as it will be the same as the end point)
	std::vector<CanvasPoint> result = {from};
	if (numberOfValues < 2) return result;
	CanvasPoint interval = {to.x - from.x, to.y - from.y, to.depth - from.depth};
	float stepX = interval.x / float(numberOfValues - 1);
	float stepY = interval.y / float(numberOfValues - 1);
	float stepD = interval.depth / float(numberOfValues - 1);

	std::vector<TexturePoint> tps = interpolateTexturePoints(from.texturePoint, to.texturePoint, numberOfValues);
	for (size_t i = 1; i < numberOfValues; i++) {
		CanvasPoint point = {from.x + float(i) * stepX, from.y + float(i) * stepY, from.depth + float(i) * stepD};
		point.texturePoint = tps[i];
		result.push_back(point);
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, size_t numberOfValues) {
	// if one or fewer values, return only the start point (as it will be the same as the end point)
	std::vector<glm::vec3> result = {from};
	if (numberOfValues < 2) return result;

	glm::vec3 interval = to - from;
	glm::vec3 step = interval / float(numberOfValues - 1);

	result.reserve(numberOfValues);
	for (size_t i = 1; i < numberOfValues; i++) {
		result.push_back(from + float(i) * step);
	}
	return result;
}


void drawLine(CanvasPoint &from, CanvasPoint &to, glm::vec3 &colour, DrawingWindow &window) {
	// CULL COMPLETELY OFF-SCREEN LINES
	if (isOffScreen(from, to)) return;

	// number of steps is one more than the maximum difference in x- or y-direction
	int numberOfValues = ceil(glm::max(abs(to.x - from.x), abs(to.y - from.y)) + 1);
	// interp along x and y, store in vector
	std::vector<CanvasPoint> line = interpolateCanvasPoints(from, to, numberOfValues);
	// pack colour
	uint32_t packedCol = (255 << 24) + (int(colour.r * 255) << 16) + (int(colour.g * 255) << 8) + int(colour.b * 255);

	// set pixels based on interpolated values (rounded)
	for (int i = 0; i < numberOfValues; i++) {
		CanvasPoint point = CanvasPoint(line[i].x, line[i].y, line[i].depth);
		// check point is on screen
		if (!isOffScreen(point)) {
			window.setPixelColour(floor(point.x), floor(point.y), packedCol);
		}
	}
}

void drawTriangle(CanvasTriangle &triangle, glm::vec3 &colour, DrawingWindow &window) {
	// CULL COMPLETELY OFF-SCREEN TRIANGLES
	if (isOffScreen(triangle.v0(), triangle.v1(), triangle.v2())) return;

	// draw lines for each side
	drawLine(triangle.v0(), triangle.v1(), colour, window);
	drawLine(triangle.v1(), triangle.v2(), colour, window);
	drawLine(triangle.v2(), triangle.v0(), colour, window);
}


CanvasPoint extraVertex(CanvasTriangle &triangle) {
	// get proportion of total height for middle vertex
	float yProp = (triangle.v1().y - triangle.v0().y) / (triangle.v2().y - triangle.v0().y);

	// get x-value and depth that are y-proportion of the way along the long edge
	float v3_x = triangle.v0().x + yProp * (triangle.v2().x - triangle.v0().x);
	float v3_y = triangle.v1().y;
	float v3_depth = triangle.v0().depth + yProp * (triangle.v2().depth - triangle.v0().depth);

	// calculate v3 texture point
	float v3_tx = triangle.v0().texturePoint.x + yProp * (triangle.v2().texturePoint.x - triangle.v0().texturePoint.x);
	float v3_ty = triangle.v0().texturePoint.y + yProp * (triangle.v2().texturePoint.y - triangle.v0().texturePoint.y);
	TexturePoint v3_texture = {v3_tx, v3_ty};

	CanvasPoint v3 = CanvasPoint(v3_x, v3_y, v3_depth);
	v3.texturePoint = v3_texture;
	return v3;
}

void fillHalfTriangle(CanvasPoint &v0, CanvasPoint &v1, CanvasPoint &v3, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window, const Material &mat = Material(), const TextureMap &texture = TextureMap()) {
	bool textured = !texture.name.empty();
	// pack colour
	uint32_t packedCol = (255 << 24) + (int(mat.diffuse.r * 255) << 16) + (int(mat.diffuse.g * 255) << 8) + int(mat.diffuse.b * 255);

	// interpolate between top/bottom and middle vertices to get left and right edges
	int minY = int(round(v0.y));
	int maxY = int(round(v1.y));
	std::vector<CanvasPoint> right = interpolateCanvasPoints(v0, v3, abs(maxY - minY) + 2); // not sure why this is +2 instead of +1
	std::vector<CanvasPoint> left = interpolateCanvasPoints(v0, v1, abs(maxY - minY) + 2);
	std::vector<CanvasPoint> row;

	// fill desired pixels
	for (int i = 0; i < left.size(); i++) {
		size_t rowLength = size_t(abs(right[i].x - left[i].x)) + 2;
		// interpolate between left and right edges to get canvas points
		row = interpolateCanvasPoints(left[i], right[i], rowLength);
		for (int j = 0; j < rowLength; j++) {
			int x = int(round(row[j].x)), y = int(round(row[j].y));
			// if within screen range and existing depth is less than new depth, draw pixel and update depth
			if (!(x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) && depthBuffer[y][x] < row[j].depth) {
				if (textured) {
					packedCol = texture.pixels[(int(round(row[j].texturePoint.y * float(texture.height))) * texture.width + int(round(row[j].texturePoint.x * float(texture.width)))) % texture.pixels.size()];
				}
				window.setPixelColour(x, y, packedCol);
				depthBuffer[y][x] = row[j].depth;
			}
		}
	}
}

void fillTriangle(CanvasTriangle &triangle, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window, const Material &mat = Material(), const TextureMap &texture = TextureMap()) {
	// CULL COMPLETELY OFF-SCREEN TRIANGLES
	//if (isOffScreen(triangle.v0(), triangle.v1(), triangle.v2())) return;

	// sort vertices by y
	for (int i = 1; i < 3; i++) {
		CanvasPoint tmpV = triangle.vertices[i];
		int j = i - 1;
		while (j >= 0 && triangle.vertices[j].y > tmpV.y) {
			triangle.vertices[j + 1] = triangle.vertices[j];
			j = j - 1;
		}
		triangle.vertices[j + 1] = tmpV;
	}

	// CALCULATE EXTRA MIDDLE VERTEX
	CanvasPoint triangle_v3 = extraVertex(triangle);

	// assign v0 through v3 based on orientation
	CanvasPoint v0 = triangle.v0(); CanvasPoint v1; CanvasPoint v2 = triangle.v2(); CanvasPoint v3;
	if (triangle_v3.x < triangle.v1().x) {v1 = triangle_v3; v3 = triangle.v1();}
	else {v1 = triangle.v1(); v3 = triangle_v3;}

	// FILL TOP HALF
	fillHalfTriangle(v0, v1, v3, depthBuffer, window, mat, texture);

	// FILL BOTTOM HALF
	fillHalfTriangle(v2, v1, v3, depthBuffer, window, mat, texture);
}


void redNoise(DrawingWindow &window) {
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			auto red = float(rand() % 256);
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void gradient(DrawingWindow &window) {
	// interp from white to black
	std::vector<float> greyVals = interpolateSingleFloats(255, 0, WIDTH);

	// set all RGB values to the stored interpolated greyValue
	for (size_t y = 0; y < HEIGHT; y++) {
		for (size_t x = 0; x < WIDTH; x++) {
			float grey = greyVals[x];
			uint32_t colour = (255 << 24) + (int(grey) << 16) + (int(grey) << 8) + int(grey);
			window.setPixelColour(x, y, colour);
		}
	}
}

void colourGradient(DrawingWindow &window) {
	// target colours
	const glm::vec3 topLeft(255, 0, 0);        // red
	const glm::vec3 topRight(0, 0, 255);       // blue
	const glm::vec3 bottomRight(0, 255, 0);    // green
	const glm::vec3 bottomLeft(255, 255, 0);   // yellow

	// interpolation vertically along left and right edge
	std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, WIDTH);
	std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);

	// for each row, interp between left and rightmost pixel
	for (size_t y = 0; y < HEIGHT; y++) {
		std::vector<glm::vec3> between = interpolateThreeElementValues(left[y], right[y], WIDTH);
		for (size_t x = 0; x < WIDTH; x++) {
			uint32_t colour = (255 << 24) + (int(between[x].r) << 16) + (int(between[x].g) << 8) + int(between[x].b);
			window.setPixelColour(x, y, colour);
		}
	}
}


void interpolationTests() {
	// TESTS
	// test single interpolation
	std::vector<float> result1;
	result1 = interpolateSingleFloats(2.2, 8.5, 7);
	for (auto i : result1) std::cout << i << " ";
	std::cout << std::endl;
	// test vec3 interpolation
	std::vector<glm::vec3> result2;
	result2 = interpolateThreeElementValues(glm::vec3(1.0, 4.0, 9.2), glm::vec3(4.0, 1.0, 9.8), 4);
	for (auto i : result2) std::cout << "(" << i.x << ", " << i.y << ", " << i.z << ") " << std::endl;
	std::cout << std::endl;
}

void textureTriangleTest(DrawingWindow &window) {
	// TEXTURE TRIANGLE TEST
	std::vector<std::vector<float>> depthBuffer(HEIGHT, std::vector<float>(WIDTH, 0));
	TextureMap texture = TextureMap("./assets/textures/triangle.ppm");
	CanvasPoint v0 = CanvasPoint(160, 10, 1);
	CanvasPoint v1 = CanvasPoint(300, 230, 1);
	CanvasPoint v2 = CanvasPoint(10, 150, 1);
	v0.texturePoint.x = 195.0f / float(texture.width); v0.texturePoint.y = 5.0f / float(texture.height);
	v1.texturePoint.x = 395.0f / float(texture.width); v1.texturePoint.y = 380.0f / float(texture.height);
	v2.texturePoint.x = 65.0f / float(texture.width); v2.texturePoint.y = 330.0f / float(texture.height);
	CanvasTriangle canvasTriangle = CanvasTriangle(v0, v1, v2);
	fillTriangle(canvasTriangle, depthBuffer, window, Material(), texture);
}


void pointCloud(std::vector<ModelTriangle> &triangles, glm::vec3 &colour, Camera &cam, DrawingWindow &window) {
	for (auto &triangle: triangles) {
		for (auto vertex : triangle.vertices) {
			// for each vertex in model triangles, project onto canvas and colour
			CanvasPoint point = projectVertexOntoCanvasPoint(vertex.position, cam);
			// check points are on-screen
			if (!isOffScreen(point)) {
				uint32_t packedCol = (255 << 24) + (int(colour.r * 255) << 16) + (int(colour.g) << 8) + int(colour.b);
				window.setPixelColour(floor(point.x), floor(point.y), packedCol);
			}
		}
	}
}

void wireFrame(std::vector<ModelTriangle> &triangles, glm::vec3 &colour, Camera &cam, DrawingWindow &window) {
	for (auto triangle: triangles) {
		// for each model triangle, project vertices onto canvas and draw resulting triangle
		CanvasPoint v0 = projectVertexOntoCanvasPoint(triangle.v0().position, cam);
		CanvasPoint v1 = projectVertexOntoCanvasPoint(triangle.v1().position, cam);
		CanvasPoint v2 = projectVertexOntoCanvasPoint(triangle.v2().position, cam);
		CanvasTriangle canvasTriangle = CanvasTriangle(v0, v1, v2);
		drawTriangle(canvasTriangle, colour, window);
	}
}

void raster(std::vector<ModelTriangle> &triangles, std::map<std::string, TextureMap> &textures, Camera &cam, DrawingWindow &window) {
	// initialise depth buffer with zero-values
	std::vector<std::vector<float>> depthBuffer(HEIGHT, std::vector<float>(WIDTH, 0));
	for (auto triangle: triangles) {
		// for each model triangle, draw a filled triangle from projected points
		CanvasPoint v0 = projectVertexOntoCanvasPoint(triangle.v0().position, cam);
		v0.texturePoint = triangle.v0().texturePoint;
		CanvasPoint v1 = projectVertexOntoCanvasPoint(triangle.v1().position, cam);
		v1.texturePoint = triangle.v1().texturePoint;
		CanvasPoint v2 = projectVertexOntoCanvasPoint(triangle.v2().position, cam);
		v2.texturePoint = triangle.v2().texturePoint;
		CanvasTriangle canvasTriangle = CanvasTriangle(v0, v1, v2);
		fillTriangle(canvasTriangle, depthBuffer, window, triangle.material, triangle.texture.empty() ? TextureMap() : textures[triangle.texture]);
	}
}

void raytrace(std::vector<ModelTriangle> &triangles, std::map<std::string, TextureMap> &textures, std::map<std::string, TextureMap> &normalMaps, std::vector<Light> &lights, Camera &cam, DrawingWindow &window) {
	omp_set_nested(1);
	#pragma omp parallel for default(shared) num_threads(6) collapse(2)
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			// get the closest intersection of a ray through current (x,y) on the image plane
			CanvasPoint canvasPoint = CanvasPoint(float(x), float(y));
			glm::vec3 rayDirection = rayFromCanvasPoint(canvasPoint, cam);
			RayTriangleIntersection surface = closestIntersection(triangles, cam.position, rayDirection, cam, MAX_RAY_BOUNCES);
			ModelTriangle &triangle = surface.intersectedTriangle;
			// if hit a triangle
			if (surface.triangleIndex != size_t(-1)) {
				glm::vec3 &colour = triangle.material.diffuse;
				float brightness = 0.0f;
				// if textured set colour from texture map
				if (!triangle.texture.empty()) {
					TextureMap &texture = textures[triangle.texture];
					TexturePoint texturePoint = barycentricTexturePoint(triangle, surface.barycentricPoint);
					uint32_t packedCol = texture.pixels[(int(round(texturePoint.y * float(texture.height))) * texture.width + int(round(texturePoint.x * float(texture.width)))) % texture.pixels.size()];
					uint32_t rgb = packedCol - (255 << 24);
					uint32_t r = rgb >> 16;
					uint32_t g = (rgb - (r << 16)) >> 8;
					uint32_t b = (rgb - (r << 16) - (g << 8));
					colour = glm::vec3(float(r) / 255.0f, float(g) / 255.0f, float(b) / 255.0f);
				}
				// if normal mapped get normal from map
				if (!triangle.normalMap.empty()) {
					TextureMap &normalMap = normalMaps[triangle.normalMap];
					TexturePoint normalPoint = barycentricNormalPoint(triangle, surface.barycentricPoint);
					uint32_t packedCol = normalMap.pixels[(int(round((normalPoint.y * float(normalMap.height))) * float(normalMap.width)) + int(round(normalPoint.x * float(normalMap.width)))) % (normalMap.width * normalMap.height)];
					uint32_t rgb = packedCol - (255 << 24);
					uint32_t r = rgb >> 16;
					uint32_t g = (rgb - (r << 16)) >> 8;
					uint32_t b = (rgb - (r << 16) - (g << 8));
					// convert from rgb [0 -> 1] to normals [-1 -> 1]
					triangle.normal = {float(r) / 255.0f * 2.0f - 1.0f, float(g) / 255.0f * 2.0f - 1.0f, float(b) / 255.0f * 2.0f - 1.0f};
				}
				// check if it is an emissive surface
				if (!triangle.material.emissive) {
					// calculate brightness (Phong, Gouraud or face-normal techniques) for each light
					for (Light light : lights) {
						// check that surface is in front of area light
						if ((light.type == AREA && dot(cross(light.e1, light.e2), surface.intersectionPoint - light.position) > 0.0f) || light.type != AREA) {
							float lightBrightness = 0.0f;
							// sample light once if point light, AREA_LIGHT_SAMPLES times if area light
							#pragma omp parallel for default(none) shared(light, triangle, triangles, surface, cam, brightnessFunction, std::cout) reduction(+:lightBrightness) num_threads(2)
							for (int s = 0; s < AREA_LIGHT_SAMPLES; s++) {
								glm::vec3 lightPoint = light.sample();
								float sampleBrightness;
								// if normal mapped, use flat shading not interpolated normals
								if (triangle.normalMap.empty()) sampleBrightness = brightnessFunction(surface, lightPoint, light.intensity, cam);
								else sampleBrightness = flatBrightness(surface, lightPoint, light.intensity, cam);
								// if in shadow set brightness to zero
								if (inShadow(triangles, surface.intersectionPoint, lightPoint, cam)) sampleBrightness = 0;
								lightBrightness += sampleBrightness;
								// exit loop after one iteration if not area light
								if (light.type != AREA) s = AREA_LIGHT_SAMPLES;
							}
							brightness += lightBrightness / float(light.type == AREA ? AREA_LIGHT_SAMPLES : 1.0f);
						}
					}
				} else brightness = 1.0f;
				// clamp between ambient lighting and one
				brightness = glm::clamp(brightness, cam.ambientLight, 1.0f);
				// adjust colour based on brightness
				colour *= brightness;
				uint32_t packedCol = (255 << 24) + (int(colour.r * 255) << 16) + (int(colour.g * 255) << 8) + int(colour.b * 255);
				window.setPixelColour(x, y, packedCol);
			}
			// otherwise set background colour
			else {
				/*
				float h = 0.5f * (rayDirection.y + 1.0f);
				glm::vec3 background = (1.0f - h) * glm::vec3(1.0, 1.0, 1.0) + h * glm::vec3(0.5, 0.7, 1.0);
				window.setPixelColour(x, y, (255 << 24) + (int(background.r * 255) << 16) + (int(background.g * 255) << 8) + int(background.b * 255));
				*/
			}
		}
	}
}


void draw(std::vector<ModelTriangle> &triangles, std::map<std::string, TextureMap> &textures, std::map<std::string, TextureMap> &normalMaps, std::vector<Light> &lights, Camera &cam, DrawingWindow &window) {
	// clear frame before drawing anything new
	window.clearPixels();

	Material red = Material({1, 0, 0});
	Material white = Material({1, 1, 1});

	Light &testLight = lights.at(0);
	// create triangle to view light source for testing
	std::vector<ModelTriangle> lightTrig;
	testLight.type == AREA ?
		lightTrig = {{
			Vertex(testLight.position),
			Vertex(testLight.e1 + testLight.position),
			Vertex(testLight.e2 + testLight.position),
			red
		}, {
			Vertex(testLight.e1 + testLight.e2 + testLight.position),
			Vertex(testLight.e1 + testLight.position),
			Vertex(testLight.e2 + testLight.position),
			red
		}}
	:
		lightTrig = {{
			Vertex(testLight.position + glm::vec3(0.0f, 0.025f, 0.0f)),
			Vertex(testLight.position + glm::vec3(-0.025f, -0.025f, 0.0f)),
			Vertex(testLight.position + glm::vec3(0.025f, -0.025f, 0.0f)),
			red
		}};
	std::vector<ModelTriangle> extras = {lightTrig};

	if (cam.orbit) {
		cam.position = cam.position * rotateY(degToRad(cam.lookSpeed));
		glm::vec3 target = {0, 0, 0};
		lookAt(target, cam);
		window.clearPixels();
	}

	switch (cam.renderMode) {
		case WIRE:
			wireFrame(triangles, white.diffuse, cam, window);
			break;
		case RASTER:
			raster(triangles, textures, cam, window);
			// render triangle where light source is
			wireFrame(extras, red.diffuse, cam, window);
			break;
		case RAYTRACE:
			raytrace(triangles, textures, normalMaps, lights, cam, window);
			break;
		default:
			break;
	}
}

void handleEvent(SDL_Event &event, std::vector<Light> &lights, Camera &cam, DrawingWindow &window) {
	Light &testLight = lights.at(0);
	if (event.type == SDL_KEYDOWN) {
		// move left
		if (event.key.keysym.sym == SDLK_a) cam.position -= normalize(cam.rotation[0]) * cam.moveSpeed;
		// move right
		else if (event.key.keysym.sym == SDLK_d) cam.position += normalize(cam.rotation[0]) * cam.moveSpeed;
		// move forwards
		else if (event.key.keysym.sym == SDLK_w) cam.position -= normalize(cam.rotation[2]) * cam.moveSpeed;
		// move backwards
		else if (event.key.keysym.sym == SDLK_s) cam.position += normalize(cam.rotation[2]) * cam.moveSpeed;
		// move up
		else if (event.key.keysym.sym == SDLK_r) cam.position += normalize(cam.rotation[1]) * cam.moveSpeed;
		// move down
		else if (event.key.keysym.sym == SDLK_f) cam.position -= normalize(cam.rotation[1]) * cam.moveSpeed;
		// orbit right
		else if (event.key.keysym.sym == SDLK_e) cam.position = cam.position * rotateY(degToRad(-cam.lookSpeed));
		// orbit left
		else if (event.key.keysym.sym == SDLK_q) cam.position = cam.position * rotateY(degToRad(cam.lookSpeed));
		// look left
		else if (event.key.keysym.sym == SDLK_LEFT) cam.rotation = rotateY(degToRad(cam.lookSpeed)) * cam.rotation;
		// look right
		else if (event.key.keysym.sym == SDLK_RIGHT) cam.rotation = rotateY(degToRad(-cam.lookSpeed)) * cam.rotation;
		// look up
		else if (event.key.keysym.sym == SDLK_UP) cam.rotation = rotateV(cam.rotation[0], degToRad(-cam.lookSpeed)) * cam.rotation;
		// look down
		else if (event.key.keysym.sym == SDLK_DOWN) cam.rotation = rotateV(cam.rotation[0], degToRad(cam.lookSpeed)) * cam.rotation;
		// reset position and rotation
		else if (event.key.keysym.sym == SDLK_SPACE) { cam.position = cam.defaultPosition; cam.rotation = cam.defaultRotation; }
		// toggle orbit
		else if (event.key.keysym.sym == SDLK_LCTRL) cam.orbit = !cam.orbit;
		// look at origin
		else if (event.key.keysym.sym == SDLK_c) { glm::vec3 target = {0, 0, 0}; lookAt(target, cam); }
		// move light forward
		else if (event.key.keysym.sym == SDLK_i) testLight.position += glm::vec3(0, 0, -1) * cam.moveSpeed;
		// move light backward
		else if (event.key.keysym.sym == SDLK_k) testLight.position += glm::vec3(0, 0, 1) * cam.moveSpeed;
		// move light left
		else if (event.key.keysym.sym == SDLK_j) testLight.position += glm::vec3(-1, 0, 0) * cam.moveSpeed;
		// move light right
		else if (event.key.keysym.sym == SDLK_l) testLight.position += glm::vec3(1, 0, 0) * cam.moveSpeed;
		// move light down
		else if (event.key.keysym.sym == SDLK_o) testLight.position += glm::vec3(0, -1, 0) * cam.moveSpeed;
		// move light up
		else if (event.key.keysym.sym == SDLK_p) testLight.position += glm::vec3(0, 1, 0) * cam.moveSpeed;
		// reset light position
		else if (event.key.keysym.sym == SDLK_u) testLight.position = testLight.defaultPosition;
		// toggle fast movement
		else if (event.key.keysym.sym == SDLK_LSHIFT) {
			float oldSpeed = cam.moveSpeed;
			cam.moveSpeed = cam.altMoveSpeed;
			cam.altMoveSpeed = oldSpeed;

			oldSpeed = cam.lookSpeed;
			cam.lookSpeed = cam.altLookSpeed;
			cam.altLookSpeed = oldSpeed;
		}
		// switch render mode
		else if (event.key.keysym.sym == SDLK_1) { cam.renderMode = WIRE; std::cout << "Wireframing..." << std::endl; }
		else if (event.key.keysym.sym == SDLK_2) { cam.renderMode = RASTER; std::cout << "Rasterising..." << std::endl; }
		else if (event.key.keysym.sym == SDLK_3) { cam.renderMode = RAYTRACE; std::cout << "Raytracing..." << std::endl; }
		// switch shading mode
		else if (event.key.keysym.sym == SDLK_4) { cam.shadingMode = FLAT; brightnessFunction = flatBrightness; std::cout << "Using flat shading..." << std::endl; }
		else if (event.key.keysym.sym == SDLK_5) { cam.shadingMode = GOURAUD; brightnessFunction = gouraudBrightness; std::cout << "Using gouraud shading..." << std::endl; }
		else if (event.key.keysym.sym == SDLK_6) { cam.shadingMode = PHONG; brightnessFunction = flatBrightness; std::cout << "Using phong shading..." << std::endl; }
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}


glm::vec3 vertexNormal(Vertex &vertex, std::vector<ModelTriangle> &triangles) {
	// for each vertex on triangles in the model, if it's the same as given vertex, add to list
	glm::vec3 total;
	for (auto &triangle : triangles) {
		for (auto &v : triangle.vertices) {
			if (v.index == vertex.index) {
				// sum the normals of neighbouring triangles
				total += triangle.normal;
			}
		}
	}

	// normalize total to get average normal
	total = normalize(total);
	return total;
}

std::tuple<std::map<std::string, Material>, std::map<std::string, TextureMap>, std::map<std::string, TextureMap>> readMtl(const std::string &filename) {
	std::map<std::string, Material> materials;
	std::map<std::string, TextureMap> textures, normals;
	std::ifstream file(filename);
	std::string line;
	std::string name;
	std::vector<std::string> splitln;
	// read from file
	while(getline(file, line)) {
		splitln = split(line, ' ');
		if (splitln[0] == "newmtl") {
			// set colour name
			name = splitln[1];
		}
		if (splitln[0] == "Kd") {
			// create new colour from given values
			float r = strtof(splitln[1].c_str(), nullptr);
			float g = strtof(splitln[2].c_str(), nullptr);
			float b = strtof(splitln[3].c_str(), nullptr);
			materials.insert({name, Material(name, glm::vec3(r, g, b))});
		}
		if (splitln[0] == "Ks") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			materials[name].mirrored = true;
		}
		if (splitln[0] == "Ke") {
			materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
			materials[name].emissive = true;
		}
		if (splitln[0] == "Ka") {
			materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
			materials[name].glassy = true;
			materials[name].refractiveIndex = strtof(splitln[1].c_str(), nullptr);
		}
		if (splitln[0] == "map_Kd") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			// find texture filename and create TextureMap object
			/*
			size_t parent = filename.find_last_of('/');
			std::string texFile = filename.substr(0, parent + 1) + splitln[1];
			*/
			std::string texFile = "./assets/textures/" + splitln[1];
			TextureMap texture = TextureMap(texFile);
			textures.insert({name, texture});
		}
		if (splitln[0] == "map_bump") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			// find texture filename and create TextureMap object
			/*size_t parent = filename.find_last_of('/');
			std::string texFile = filename.substr(0, parent + 1) + splitln[1];
			*/
			std::string texFile = "./assets/normals/" + splitln[1];
			TextureMap normal = TextureMap(texFile);
			normals.insert({name, normal});
		}
	}

	return {materials, textures, normals};
}

std::tuple<std::vector<ModelTriangle>, std::map<std::string, TextureMap>, std::map<std::string, TextureMap>> readObj(const std::string &filename, float modelScale = 0.35) {
	std::vector<ModelTriangle> triangles;
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> texturePoints;
	std::map<std::string, Material> materials;
	std::map<std::string, TextureMap> textures, normalMaps;
	TextureMap texture, normalMap;
	Material material;
	std::string line;
	std::vector<std::string> splitln;
	std::ifstream file(filename);
	// read from file
	while(getline(file, line)) {
		// split line into segments
		splitln = split(line, ' ');
		// read from mtl file
		if (splitln[0] == "mtllib") {
			size_t parent = filename.find_last_of('/');
			std::string matFile = filename.substr(0, parent + 1) + splitln[1];
			std::tuple<std::map<std::string, Material>, std::map<std::string, TextureMap>, std::map<std::string, TextureMap>> mtl = readMtl(matFile);
			materials = std::get<0>(mtl);
			textures = std::get<1>(mtl);
			normalMaps = std::get<2>(mtl);
		}
		if (splitln[0] == "usemtl") {
			// look up colour from palette
			material = materials[splitln[1]];
			// check if texture and normal map exist for this material, otherwise reset them
			auto it = textures.find(splitln[1]);
			if (it != textures.end()) texture = it->second;
			else texture = TextureMap();
			it = normalMaps.find(splitln[1]);
			if (it != normalMaps.end()) normalMap = it->second;
			else normalMap = TextureMap();
		}
		if (splitln[0] == "v") {
			// convert coords to floats and scale by 'modelScale' (ignoring first index which will be 'v')
			glm::vec3 vertex = modelScale * glm::vec3(strtod(splitln[1].c_str(), nullptr), strtod(splitln[2].c_str(), nullptr), strtod(splitln[3].c_str(), nullptr));
			vertices.push_back(vertex);
		}
		if (splitln[0] == "vt") {
			// convert to floats and create texture point (ignoring first index which will be 'vt')
			TexturePoint texturePoint = {float(strtod(splitln[1].c_str(), nullptr)), float(strtod(splitln[2].c_str(), nullptr))};
			texturePoints.push_back(texturePoint);
		}
		if (splitln[0] == "f") {
			// split line into indices
			std::vector<std::string> vIndices, tIndices;
			for (size_t i = 1; i < splitln.size(); i++) {
				std::vector<std::string> vt = split(splitln[i], '/');
				vIndices.push_back(vt[0]);
				if (vt.size() > 1) tIndices.push_back(vt[1]);
			}
			// convert indices to ints and zero-index (ignoring first char 'f')
			int iv0 = int(strtol(vIndices[0].c_str(), nullptr, 10) - 1);
			int iv1 = int(strtol(vIndices[1].c_str(), nullptr, 10) - 1);
			int iv2 = int(strtol(vIndices[2].c_str(), nullptr, 10) - 1);
			// set vertices and their indices
			Vertex v0 = Vertex(vertices[iv0]), v1 = Vertex(vertices[iv1]), v2 = Vertex(vertices[iv2]);
			v0.index = iv0; v1.index = iv1; v2.index = iv2;
			// set texture points if used
			if (!tIndices.empty()) {
				int it0 = int(strtol(tIndices[0].c_str(), nullptr, 10) - 1);
				int it1 = int(strtol(tIndices[1].c_str(), nullptr, 10) - 1);
				int it2 = int(strtol(tIndices[2].c_str(), nullptr, 10) - 1);
				TexturePoint t0 = texturePoints[it0], t1 = texturePoints[it1], t2 = texturePoints[it2];
				v0.texturePoint = t0; v1.texturePoint = t1; v2.texturePoint = t2;
			}
			// set normal map points if used
			if (!normalMap.name.empty()) {
				v0.normalPoint = TexturePoint(0, 0);
				v1.normalPoint = TexturePoint((v1.position - v0.position).x, (v1.position - v0.position).y);
				v2.normalPoint = TexturePoint((v2.position - v0.position).x, (v2.position - v0.position).y);
			}
			// calculate triangle face normal
			glm::vec3 normal = normalize(cross(v1.position - v0.position, v2.position - v0.position));
			// create new triangle
			ModelTriangle modelTriangle = {v0, v1, v2, material};
			modelTriangle.normal = normal;
			modelTriangle.texture = texture.name.empty() ? "" : material.name;
			modelTriangle.normalMap = normalMap.name.empty() ? "" : material.name;

			// add triangle
			triangles.emplace_back(modelTriangle);
		}
	}

	// calculate vertex normals based on neighbouring triangles
	for (auto &triangle : triangles) {
		for (auto &v : triangle.vertices) {
			v.normal = vertexNormal(v, triangles);
		}
	}

	return {triangles, textures, normalMaps};
}


void animate(DrawingWindow &window) {
	// READ OBJ AND MTL FILES
	// create vector of triangles and textures
	std::vector<ModelTriangle> triangles;
	std::map<std::string, TextureMap> textures, normalMaps;

	// read obj and mtl files and store
	auto scene = readObj(SCENE_TWO);
	triangles = std::get<0>(scene);
	textures = std::get<1>(scene);
	normalMaps = std::get<2>(scene);

	// INITIALISE LIGHTS AND CAMERA
	Camera cam = Camera(glm::vec3(0, 0, 3), 2.0f);
	Light ceiling = Light(AREA, 20.0f);
	Light pointLight = Light(POINT, glm::vec3(0.0, 0.75, 0.5), 10.0f);
	std::vector<Light> lights = {pointLight};
	cam.shadingMode = FLAT;
	brightnessFunction = flatBrightness;

	glm::vec3 white = {1, 1, 1};

	system("rm -r ./animation");
	system("mkdir ./animation");
	int frames = 0, sectionOne = 120, sectionTwo = 180, sectionThree = 360, sectionFour = 480, sectionFive = 600;

	for (int frame = frames; frame < sectionOne; frame++) {
		frames = frame + 1;
		window.clearPixels();
		wireFrame(triangles, white, cam, window);
		window.renderFrame();

		cam.position = cam.position * rotateY(degToRad(360.0f / float(sectionOne)));
		glm::vec3 target = {0, 0, 0};
		lookAt(target, cam);

		std::stringstream filename;
		filename << "./animation/" << std::to_string(frame) << ".ppm";
		window.savePPM(filename.str());
	}

	for (int frame = frames; frame < sectionTwo; frame++) {
		frames = frame + 1;
		window.clearPixels();
		raster(triangles, textures, cam, window);
		window.renderFrame();

		if (frame < sectionOne + (sectionTwo - sectionOne) / 2) {
			cam.position += cam.moveSpeed * glm::vec3(0, 0, 2);
			cam.rotation = rotateV(glm::vec3(0, 0, -1), degToRad(cam.altLookSpeed)) * cam.rotation;
		} else {
			cam.position += cam.moveSpeed * glm::vec3(0, 0, -2);
			cam.rotation = rotateV(glm::vec3(0, 0, -1), degToRad(cam.altLookSpeed)) * cam.rotation;
		}

		glm::vec3 target = {0, 0, 0};
		lookAt(target, cam);

		std::stringstream filename;
		filename << "./animation/" << std::to_string(frame) << ".ppm";
		window.savePPM(filename.str());
	}

	// disable ceiling light brightness
	triangles.at(8).material.emissive = false;
	triangles.at(9).material.emissive = false;
	lights.at(0).intensity = 0.0f;

	for (int frame = frames; frame < sectionThree; frame++) {
		frames = frame + 1;
		window.clearPixels();
		raytrace(triangles, textures, normalMaps, lights, cam, window);
		window.renderFrame();

		if (frame < sectionThree + (sectionThree - sectionTwo) / 2) {
			lights.at(0).position = lights.at(0).position * rotateY(degToRad(360.0f / float(sectionOne)));
			lights.at(0).intensity += 10.0f / (float(sectionThree - sectionTwo) / 2.0f);
		} else {
			lights.at(0).position += 20 / float(sectionThree - sectionTwo) * glm::vec3(-0.1, 0, -0.1);
			lights.at(0).intensity -= 10.0f / (float(sectionThree - sectionTwo) / 2.0f);
		}

		std::stringstream filename;
		filename << "./animation/" << std::to_string(frame) << ".ppm";
		window.savePPM(filename.str());
	}

	// enable ceiling light brightness
	triangles.at(8).material.emissive = true;
	triangles.at(9).material.emissive = true;
	lights.at(0) = ceiling;
	lights.at(0).intensity = 0.0f;
	cam.shadingMode = GOURAUD;
	brightnessFunction = gouraudBrightness;

	for (int frame = frames; frame < sectionFour; frame++) {
		frames = frame + 1;
		window.clearPixels();
		raytrace(triangles, textures, normalMaps, lights, cam, window);
		window.renderFrame();

		if (frame < sectionThree + 30) {
			lights.at(0).intensity += 20.0f / 30.0f;
		}

		std::stringstream filename;
		filename << "./animation/" << std::to_string(frame) << ".ppm";
		window.savePPM(filename.str());
	}

	cam.shadingMode = PHONG;
	brightnessFunction = flatBrightness;

	scene = readObj(SCENE_ONE);
	triangles = std::get<0>(scene);
	textures = std::get<1>(scene);
	normalMaps = std::get<2>(scene);

	for (int frame = frames; frame < sectionFive; frame++) {
		frames = frame + 1;
		window.clearPixels();
		raytrace(triangles, textures, normalMaps, lights, cam, window);
		window.renderFrame();

		std::stringstream filename;
		filename << "./animation/" << std::to_string(frame) << ".ppm";
		window.savePPM(filename.str());
	}
}


[[noreturn]] int main(int argc, char *argv[]) {
	// seed rand
	srand(time(nullptr));

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	// READ OBJ AND MTL FILES
	// create vector of triangles and textures
	std::vector<ModelTriangle> triangles;
	std::map<std::string, TextureMap> textures, normalMaps;

	// read obj and mtl files and store
	auto scene = readObj(SCENE_TWO);
	triangles.insert(triangles.end(), std::get<0>(scene).begin(), std::get<0>(scene).end());
	textures.insert(std::get<1>(scene).begin(), std::get<1>(scene).end());
	normalMaps.insert(std::get<2>(scene).begin(), std::get<2>(scene).end());

	// INITIALISE LIGHTS AND CAMERA
	std::vector<Light> lights;

	Camera camera = Camera(glm::vec3(0, 0, 3), 2.0f);
	Light ceiling = Light(AREA, 20.0f);
	Light pointLight = Light(POINT, glm::vec3(0.0, 0.5, 0.5), 10.0f);

	lights.insert(lights.end(), {pointLight});
	brightnessFunction = flatBrightness;

	// UNCOMMENT FOR ANIMATION
	//animate(window);

	time_t start = time(nullptr);
	time_t elapsed;
	int frames = 0;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, lights, camera, window);
		draw(triangles, textures, normalMaps, lights, camera, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();

		frames++;
		elapsed = time(nullptr) - start;
		// basic frame time calculator
		if (elapsed >= 1) {
			float ms = float(elapsed) * 1000.0f / float(frames);
			if (ms >= 60000.0f) { ms /= 1000.0f; std::cout << "frame time: " << int(ms / 60) << "m " << int(ms) % 60 << "s" << std::endl; }
			else if (ms >= 1000.0f) { ms /= 1000.0f; std::cout << "frame time: " << ms << "s" << std::endl; }
			else std::cout << "frame time: " << ms << "ms" << std::endl;

			start = time(nullptr);
			frames = 0;
		}
	}
}