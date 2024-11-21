#include <algorithm>
#include <map>
#include <vector>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <TextureMap.h>
#include <glm/glm.hpp>
#include <Camera.h>
#include <PointLight.h>
#include <Light.h>
#include <QuadLight.h>

#define WIDTH 640
#define HEIGHT 480

#define DEG_RAD (2.0f * M_PI / 180.0f)
#define MIN_INTERSECT_DISTANCE 0.001f
#define RENDER_TYPES 3
#define AMBIENT 0.2f

const std::string BOX_OBJ = "./assets/cornell-box/cornell-box.obj";
const std::string SPH_OBJ = "./assets/sphere/sphere.obj";
const std::string HACK_OBJ = "./assets/hackspace-logo/logo.obj";

bool ORBIT = false;
int RENDER = 0;

RayTriangleIntersection closestIntersection(std::vector<ModelTriangle> &triangles, glm::vec3 &startPoint, glm::vec3 &rayDirection) {
	// normalise ray direction and initialise closest distance
	rayDirection = normalize(rayDirection);
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
		if (t > MIN_INTERSECT_DISTANCE && 1 / t > inverseClosestDistance && u >= 0 && u <= 1.0 && v >= 0 && v <= 1.0 && (u + v) <= 1.0) {
			intersection = RayTriangleIntersection(startPoint + t * rayDirection, t, triangle, i, {u, v, 1.0f - u - v});
			inverseClosestDistance = 1 / t;
		}
	}

	if (intersection.intersectedTriangle.colour.mirrored) {
		glm::vec3 incidence = rayDirection;
		glm::vec3 normal = intersection.intersectedTriangle.normal;
		glm::vec3 start = intersection.intersectionPoint;
		glm::vec3 reflect = incidence - 2.0f * normal * dot(incidence, normal);
		return closestIntersection(triangles, start, reflect);
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

bool inShadow(std::vector<ModelTriangle> &triangles, glm::vec3 &surface, Light &light) {
	// calculate shadowRay direction
	glm::vec3 shadowRay = light.position - surface;
	// get the closest intersection of shadowRay from surface
	RayTriangleIntersection obstacle = closestIntersection(triangles, surface, shadowRay);
	// if obstacle doesn't exist or is further than light source - no shadow
	if (obstacle.triangleIndex == size_t(-1)) return false;
	if (distance(surface, light.position) < distance(surface, obstacle.intersectionPoint)) return false;
	return true;
}

float proximityLighting(glm::vec3 &point, Light &light) {
	float dist = distance(point, light.position);
	// + 1 in denominator to avoid infinite brightness at light source
	float brightness = light.intensity / float(4 * M_PI * dist * dist + 1);
	return brightness;
}

float angularLighting(glm::vec3 &point, glm::vec3 &normal, Light &light) {
	float angle = dot(normal, normalize(light.position - point));
	float brightness = angle;
	return brightness;
}

float specularLighting(glm::vec3 &point, glm::vec3 &normal, Light &light, glm::vec3 &camPos, float n = 256) {
	glm::vec3 view = normalize(camPos - point);
	glm::vec3 incidence = normalize(point - light.position);
	glm::vec3 reflection = incidence - 2.0f * normal * dot(incidence, normal);
	float brightness = pow(dot(reflection, view), n);
	brightness = glm::max(brightness, 0.0f);
	return brightness;
}

float surfaceBrightness(glm::vec3 &point, glm::vec3 &normal, Light &light, Camera &cam) {
	// calculate proximity lighting
	float prox = proximityLighting(point, light);
	// calculate angle of incidence lighting
	float ang = angularLighting(point, normal, light);
	// calculate specular lighting
	float spec = specularLighting(point, normal, light, cam.position);

	float brightness = prox * ang + spec;
	brightness = glm::clamp(brightness, 0.0f, 1.0f);
	return brightness;
}

float gouraudBrightness(RayTriangleIntersection &intersection, Light &light, Camera &cam) {
	glm::vec3 &bCoords = intersection.barycentricPoint;
	glm::vec3 &v0_norm = intersection.intersectedTriangle.v0().normal;
	glm::vec3 &v1_norm = intersection.intersectedTriangle.v1().normal;
	glm::vec3 &v2_norm = intersection.intersectedTriangle.v2().normal;
	float v0_brightness = surfaceBrightness(intersection.intersectionPoint, v0_norm, light, cam);
	float v1_brightness = surfaceBrightness(intersection.intersectionPoint, v1_norm, light, cam);
	float v2_brightness = surfaceBrightness(intersection.intersectionPoint, v2_norm, light, cam);
	float brightness = bCoords.x * v1_brightness + bCoords.y * v2_brightness + bCoords.z * v0_brightness;
	return brightness;
}

float phongBrightness(RayTriangleIntersection &intersection, Light &light, Camera &cam) {
	ModelTriangle &triangle = intersection.intersectedTriangle;
	glm::vec3 &bCoords = intersection.barycentricPoint;
	glm::vec3 normal = normalize(bCoords.x * triangle.v1().normal + bCoords.y * triangle.v2().normal + bCoords.z * triangle.v0().normal);
	float brightness = surfaceBrightness(intersection.intersectionPoint, normal, light, cam);
	return brightness;
}


TexturePoint interpolateTexturePoint(ModelTriangle &triangle, glm::vec3 &bCoords) {
	TexturePoint tp = bCoords.x * triangle.v1().texturePoint + bCoords.y * triangle.v2().texturePoint + bCoords.z * triangle.v0().texturePoint;
	return tp;
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
	if (numberOfValues < 2) return std::vector<float>({from});
	float interval = to - from;
	float step = interval / float(numberOfValues - 1);

	std::vector<float> result;
	result.reserve(numberOfValues);
	for (size_t i = 0; i < numberOfValues; i++) {
		result.push_back(from + float(i) * step);
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, size_t numberOfValues) {
	// if one or fewer values, return only the start point (as it will be the same as the end point)
	if (numberOfValues < 2) return std::vector<glm::vec3>({from});
	glm::vec3 interval = to - from;
	glm::vec3 step = interval / float(numberOfValues - 1);

	std::vector<glm::vec3> result;
	result.reserve(numberOfValues);
	for (size_t i = 0; i < numberOfValues; i++) {
		result.push_back(from + float(i) * step);
	}
	return result;
}


void drawLine(CanvasPoint &from, CanvasPoint &to, Colour &colour, DrawingWindow &window) {
	// CULL COMPLETELY OFF-SCREEN LINES
	if (isOffScreen(from, to)) return;

	// number of steps is one more than the maximum difference in x- or y-direction
	int numberOfValues = ceil(glm::max(abs(to.x - from.x), abs(to.y - from.y)) + 1);
	// interp along x and y, store in vector
	std::vector<glm::vec3> line = interpolateThreeElementValues(glm::vec3(from.x, from.y, from.depth), glm::vec3(to.x, to.y, to.depth), numberOfValues);
	// pack colour
	uint32_t packedCol = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);

	// set pixels based on interpolated values (rounded)
	for (size_t i = 0; i < numberOfValues; i++) {
		CanvasPoint point = CanvasPoint(line[i].x, line[i].y, line[i].z);
		// check point is on screen and not behind camera
		if (!isOffScreen(point) && point.depth > 0) {
			window.setPixelColour(floor(point.x), floor(point.y), packedCol);
		}
	}
}

void drawTriangle(CanvasTriangle &triangle, Colour &colour, DrawingWindow &window) {
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
	// interpolate Xs along desired edge and get x value that is 'y proportion' of the way along
	std::vector<float> edgeXs = interpolateSingleFloats(triangle.v0().x, triangle.v2().x, size_t(triangle.v2().y - triangle.v0().y + 1));
	float v3_x = edgeXs[ceil(yProp * float(edgeXs.size() - 1))];
	// calculate v3 depth
	std::vector<float> edgeDepth = interpolateSingleFloats(triangle.v0().depth, triangle.v2().depth, size_t(triangle.v2().y - triangle.v0().y + 1));
	float v3_depth = edgeDepth[ceil(yProp * float(edgeDepth.size() - 1))];
	// calculate v3 texture point
	std::vector<glm::vec3> edgeTexture = interpolateThreeElementValues(
		glm::vec3(triangle.v0().texturePoint.x, triangle.v0().texturePoint.y, 0),
		glm::vec3(triangle.v2().texturePoint.x, triangle.v2().texturePoint.y, 0),
		size_t(triangle.v2().y - triangle.v0().y + 1));
	glm::vec3 v3_texture = edgeTexture[ceil(yProp * float(edgeDepth.size() - 1))];

	CanvasPoint v3 = CanvasPoint(v3_x, triangle.v1().y, v3_depth);
	v3.texturePoint = TexturePoint(v3_texture.x, v3_texture.y);
	return v3;
}

void fillHalfTriangle(bool top, CanvasPoint &v0, CanvasPoint &v1, CanvasPoint &v2, CanvasPoint &v3, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window, const Colour &colour = Colour(), const TextureMap &texture = TextureMap()) {
	bool textured = !texture.name.empty();
	int minY, minX, maxY, maxX;
	std::vector<float> upperBound, lowerBound, rightDepth, leftDepth, rowDepth;
	std::vector<glm::vec3> rightTexture, leftTexture, rowTexture;
	// pack colour
	uint32_t packedCol = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;

	if (top) {
		// interpolate along top edges to get x-boundaries and depth values
		minY = floor(v0.y);
		maxY = ceil(v1.y);
		minX = floor(glm::min(v1.x, v0.x));
		maxX = ceil(glm::max(v3.x, v0.x));
		upperBound = interpolateSingleFloats(v0.x, v3.x, maxY - minY + 1);
		lowerBound = interpolateSingleFloats(v0.x, v1.x, maxY - minY + 1);
		rightDepth = interpolateSingleFloats(v0.depth, v3.depth, maxY - minY + 1);
		leftDepth = interpolateSingleFloats(v0.depth, v1.depth, maxY - minY + 1);
		if (textured) {
			rightTexture = interpolateThreeElementValues(
				glm::vec3(v0.texturePoint.x, v0.texturePoint.y, 0),
				glm::vec3(v3.texturePoint.x, v3.texturePoint.y, 0),
				maxY - minY + 1);
			leftTexture = interpolateThreeElementValues(
				glm::vec3(v0.texturePoint.x, v0.texturePoint.y, 0),
				glm::vec3(v1.texturePoint.x, v1.texturePoint.y, 0),
				maxY - minY + 1);
		}
	} else {
		// interpolate along bottom edges to get x-boundaries and depth values
		minY = floor(v1.y);
		maxY = ceil(v2.y);
		minX = floor(glm::min(v1.x, v2.x));
		maxX = ceil(glm::max(v3.x, v2.x));
		upperBound = interpolateSingleFloats(v3.x, v2.x, maxY - minY + 1);
		lowerBound = interpolateSingleFloats(v1.x, v2.x, maxY - minY + 1);
		rightDepth = interpolateSingleFloats(v3.depth, v2.depth, maxY - minY + 1);
		leftDepth = interpolateSingleFloats(v1.depth, v2.depth, maxY - minY + 1);
		if (textured) {
			rightTexture = interpolateThreeElementValues(
				glm::vec3(v3.texturePoint.x, v3.texturePoint.y, 0),
				glm::vec3(v2.texturePoint.x, v2.texturePoint.y, 0),
				maxY - minY + 1);
			leftTexture = interpolateThreeElementValues(
				glm::vec3(v1.texturePoint.x, v1.texturePoint.y, 0),
				glm::vec3(v2.texturePoint.x, v2.texturePoint.y, 0),
				maxY - minY + 1);
		}
	}

	// fill desired pixels
	for (int y = minY; y < maxY; y++) {
		// check if within screen height
		if (!(y < 0 || y >= HEIGHT)) {
			// interpolate between left and right edges to get pixel depths
			rowDepth = interpolateSingleFloats(leftDepth[y - minY], rightDepth[y - minY], ceil(upperBound[y - minY] - lowerBound[y - minY] + 1));
			if (textured) rowTexture = interpolateThreeElementValues(leftTexture[y - minY], rightTexture[y - minY], ceil(upperBound[y - minY] - lowerBound[y - minY] + 1));
			for (int x = minX; x < maxX; x++) {
				// if within triangle range and screen range
				if (float(x) >= lowerBound[y - minY] && float(x) <= upperBound[y - minY] && !(x < 0 || x >= WIDTH)) {
					// if existing depth is less than current depth, draw pixel and update depth
					if (depthBuffer.at(y).at(x) < rowDepth[x - int(ceil(lowerBound[y - minY]))]) {
						if (textured) {
							glm::vec3 texturePoint = rowTexture[x - int(ceil(lowerBound[y - minY]))];
							packedCol = texture.pixels[int(round(texturePoint.y * float(texture.height))) * texture.width + int(round(texturePoint.x * float(texture.width)))];
						}
						window.setPixelColour(x, y, packedCol);
						depthBuffer.at(y).at(x) = rowDepth[x - int(ceil(lowerBound[y - minY]))];
					}
				}
			}
		}
	}
}

void fillTriangle(CanvasTriangle &triangle, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window, const Colour &colour = Colour(), const TextureMap &textureMap = TextureMap()) {
	// CULL COMPLETELY OFF-SCREEN TRIANGLES
	if (isOffScreen(triangle.v0(), triangle.v1(), triangle.v2())) return;

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
	fillHalfTriangle(true, v0, v1, v2, v3, depthBuffer, window, colour, textureMap);

	// FILL BOTTOM HALF
	fillHalfTriangle(false, v0, v1, v2, v3, depthBuffer, window, colour, textureMap);
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
	TextureMap texture = TextureMap("./assets/triangle/texture.ppm");
	CanvasPoint v0 = CanvasPoint(160, 10, 1);
	CanvasPoint v1 = CanvasPoint(300, 230, 1);
	CanvasPoint v2 = CanvasPoint(10, 150, 1);
	v0.texturePoint.x = 195.0f / float(texture.width); v0.texturePoint.y = 5.0f / float(texture.height);
	v1.texturePoint.x = 395.0f / float(texture.width); v1.texturePoint.y = 380.0f / float(texture.height);
	v2.texturePoint.x = 65.0f / float(texture.width); v2.texturePoint.y = 330.0f / float(texture.height);
	CanvasTriangle canvasTriangle = CanvasTriangle(v0, v1, v2);
	fillTriangle(canvasTriangle, depthBuffer, window, Colour(), texture);
}


void pointCloud(std::vector<ModelTriangle> &triangles, Camera &cam, DrawingWindow &window) {
	for (auto &triangle: triangles) {
		for (auto vertex : triangle.vertices) {
			// for each vertex in model triangles, project onto canvas and colour white
			CanvasPoint point = projectVertexOntoCanvasPoint(vertex.position, cam);
			// check points are on-screen
			if (!isOffScreen(point)) {
				uint32_t packedCol = (255 << 24) + (255 << 16) + (255 << 8) + 255;
				window.setPixelColour(size_t(point.x), size_t(point.y), packedCol);
			}
		}
	}
}

void wireFrame(std::vector<ModelTriangle> &triangles, Camera &cam, DrawingWindow &window) {
	Colour colour = Colour(255, 255, 255);
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
		fillTriangle(canvasTriangle, depthBuffer, window, triangle.colour, textures[triangle.texture]);
	}
}

void raytrace(std::vector<ModelTriangle> &triangles, std::map<std::string, TextureMap> &textures, std::vector<std::reference_wrapper<Light>> &lights, Camera &cam, DrawingWindow &window) {
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			// get the closest intersection of a ray through current (x,y) on the image plane
			CanvasPoint canvasPoint = CanvasPoint(float(x), float(y));
			glm::vec3 rayDirection = rayFromCanvasPoint(canvasPoint, cam);
			RayTriangleIntersection surface = closestIntersection(triangles, cam.position, rayDirection);
			if (surface.triangleIndex != size_t(-1)) {
				Colour colour = surface.intersectedTriangle.colour;
				std::cout << (colour.mirrored ? "yes" : "");
				float brightness = 0.0f;
				// if textured get colour from texture map
				if (!surface.intersectedTriangle.texture.empty()) {
					TextureMap &texture = textures[surface.intersectedTriangle.texture];
					TexturePoint texturePoint = interpolateTexturePoint(surface.intersectedTriangle, surface.barycentricPoint);
					uint32_t packedCol = texture.pixels[int(round(texturePoint.y * float(texture.height))) * texture.width + int(round(texturePoint.x * float(texture.width)))];
					uint32_t rgb = packedCol - (255 << 24);
					uint32_t r = rgb >> 16;
					uint32_t g = (rgb - (r << 16)) >> 8;
					uint32_t b = (rgb - (r << 16) - (g << 8));
					colour = Colour(int(r), int(g), int(b));
				}
				// calculate brightness (Phong, Gouraud or face-normal techniques) for each light
				for (Light light : lights) {
					float newBrightness = gouraudBrightness(surface, light, cam);
					// if in shadow set brightness to zero
					if (inShadow(triangles, surface.intersectionPoint, light)) newBrightness = 0;
					brightness += newBrightness;
				}
				// clamp between ambient lighting and one
				brightness = glm::clamp(brightness, AMBIENT, 1.0f);
				// adjust colour based on brightness
				colour = {
					int(float(colour.red) * brightness),
					int(float(colour.green) * brightness),
					int(float(colour.blue) * brightness)
				};
				uint32_t packedCol = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
				window.setPixelColour(x, y, packedCol);
			}
		}
	}
}


void draw(std::vector<ModelTriangle> &triangles, std::map<std::string, TextureMap> &textures, std::vector<std::reference_wrapper<Light>> &lights, Camera &cam, DrawingWindow &window) {
	Light &testLight = lights.at(0);
	// create triangle for light source
	ModelTriangle lightTrig = {
		{testLight.position.x - 0.025f, testLight.position.y - 0.025f, testLight.position.z},
		{testLight.position.x + 0.025f, testLight.position.y - 0.025f, testLight.position.z},
		{testLight.position.x, testLight.position.y + 0.025f, testLight.position.z},
		Colour(255, 255, 255)
	};
	std::vector<ModelTriangle> extras = {lightTrig};

	if (ORBIT) {
		cam.position = cam.position * rotateY(DEG_RAD * cam.lookSpeed);
		glm::vec3 target = {0, 0, 0};
		lookAt(target, cam);
		window.clearPixels();
	}

	switch (RENDER) {
		case 0:
			window.clearPixels();
			wireFrame(triangles, cam, window);
			break;
		case 1:
			window.clearPixels();
			raster(triangles, textures, cam, window);
			// render triangle where light source is
			wireFrame(extras, cam, window);
			break;
		case 2:
			window.clearPixels();
			raytrace(triangles, textures, lights, cam, window);
			// render triangle where light source is
			wireFrame(extras, cam, window);
			break;
		default:
			break;
	}
}

void handleEvent(SDL_Event &event, std::vector<std::reference_wrapper<Light>> &lights, Camera &cam, DrawingWindow &window) {
	Light &testLight = lights.at(0);
	if (event.type == SDL_KEYDOWN) {
		// move left
		if (event.key.keysym.sym == SDLK_a) {
			cam.position -= normalize(cam.rotation[0]) * cam.moveSpeed;
			window.clearPixels();
		}
		// move right
		else if (event.key.keysym.sym == SDLK_d) {
			cam.position += normalize(cam.rotation[0]) * cam.moveSpeed;
			window.clearPixels();
		}
		// move forwards
		else if (event.key.keysym.sym == SDLK_w) {
			cam.position -= normalize(cam.rotation[2]) * cam.moveSpeed;
			window.clearPixels();
		}
		// move backwards
		else if (event.key.keysym.sym == SDLK_s) {
			cam.position += normalize(cam.rotation[2]) * cam.moveSpeed;
			window.clearPixels();
		}
		// move up
		else if (event.key.keysym.sym == SDLK_r) {
			cam.position += normalize(cam.rotation[1]) * cam.moveSpeed;
			window.clearPixels();
		}
		// move down
		else if (event.key.keysym.sym == SDLK_f) {
			cam.position -= normalize(cam.rotation[1]) * cam.moveSpeed;
			window.clearPixels();
		}
		// orbit right
		else if (event.key.keysym.sym == SDLK_e) {
			cam.position = cam.position * rotateY(DEG_RAD * -cam.lookSpeed);
			window.clearPixels();
		}
		// orbit left
		else if (event.key.keysym.sym == SDLK_q) {
			cam.position = cam.position * rotateY(DEG_RAD * cam.lookSpeed);
			window.clearPixels();
		}
		// look left
		else if (event.key.keysym.sym == SDLK_LEFT) {
			cam.rotation = rotateY(DEG_RAD * cam.lookSpeed) * cam.rotation;
			window.clearPixels();
		}
		// look right
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			cam.rotation = rotateY(DEG_RAD * -cam.lookSpeed) * cam.rotation;
			window.clearPixels();
		}
		// look up
		else if (event.key.keysym.sym == SDLK_UP) {
			cam.rotation = rotateV(cam.rotation[0], DEG_RAD * -cam.lookSpeed) * cam.rotation;
			window.clearPixels();
		}
		// look down
		else if (event.key.keysym.sym == SDLK_DOWN) {
			cam.rotation = rotateV(cam.rotation[0], DEG_RAD * cam.lookSpeed) * cam.rotation;
			window.clearPixels();
		}
		// reset position and rotation
		else if (event.key.keysym.sym == SDLK_SPACE) {
			cam.position = cam.defaultPosition;
			cam.rotation = cam.defaultRotation;
			window.clearPixels();
		}
		// toggle orbit
		else if (event.key.keysym.sym == SDLK_g) {
			ORBIT = !ORBIT;
		}
		// look at origin
		else if (event.key.keysym.sym == SDLK_t) {
			glm::vec3 target = {0, 0, 0};
			lookAt(target, cam);
			window.clearPixels();
		}
		// move light forward
		else if (event.key.keysym.sym == SDLK_i) {
			testLight.position += glm::vec3(0, 0, -1) * cam.moveSpeed;
		}
		// move light backward
		else if (event.key.keysym.sym == SDLK_k) {
			testLight.position += glm::vec3(0, 0, 1) * cam.moveSpeed;
		}
		// move light left
		else if (event.key.keysym.sym == SDLK_j) {
			testLight.position += glm::vec3(-1, 0, 0) * cam.moveSpeed;
		}
		// move light right
		else if (event.key.keysym.sym == SDLK_l) {
			testLight.position += glm::vec3(1, 0, 0) * cam.moveSpeed;
		}
		// move light down
		else if (event.key.keysym.sym == SDLK_o) {
			testLight.position += glm::vec3(0, -1, 0) * cam.moveSpeed;
		}
		// move light up
		else if (event.key.keysym.sym == SDLK_p) {
			testLight.position += glm::vec3(0, 1, 0) * cam.moveSpeed;
		}
		// reset light position
		else if (event.key.keysym.sym == SDLK_u) {
			testLight.position = testLight.defaultPosition;
		}
		// toggle fast movement
		else if (event.key.keysym.sym == SDLK_LSHIFT) {
			float oldSpeed = cam.moveSpeed;
			cam.moveSpeed = cam.altMoveSpeed;
			cam.altMoveSpeed = oldSpeed;

			oldSpeed = cam.lookSpeed;
			cam.lookSpeed = cam.altLookSpeed;
			cam.altLookSpeed = oldSpeed;
		}
		// rotate through render modes
		else if (event.key.keysym.sym == SDLK_LCTRL) {
			RENDER++;
			RENDER %= (RENDER_TYPES);
		}
		// unfilled triangle
		else if (event.key.keysym.sym == SDLK_h) {
			Colour colour = Colour(rand() % 255, rand() % 255, rand() % 255);
			CanvasTriangle canvasTriangle = CanvasTriangle(
					CanvasPoint(float(rand() % WIDTH), float(rand() % HEIGHT), 1),
					CanvasPoint(float(rand() % WIDTH), float(rand() % HEIGHT), 1),
					CanvasPoint(float(rand() % WIDTH), float(rand() % HEIGHT), 1));
			drawTriangle(canvasTriangle, colour, window);
		}
		// filled triangle
		else if (event.key.keysym.sym == SDLK_y) {
			std::vector<std::vector<float>> depthBuffer(HEIGHT, std::vector<float>(WIDTH, 0));
			Colour colour = Colour(rand() % 255, rand() % 255, rand() % 255);
			Colour white = Colour(255, 255, 255);
			CanvasTriangle canvasTriangle = CanvasTriangle(
					CanvasPoint(float(rand() % WIDTH), float(rand() % HEIGHT), 1),
					CanvasPoint(float(rand() % WIDTH), float(rand() % HEIGHT), 1),
					CanvasPoint(float(rand() % WIDTH), float(rand() % HEIGHT), 1));
			fillTriangle(canvasTriangle, depthBuffer, window, colour);
			drawTriangle(canvasTriangle, white, window);
		}
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


std::pair<std::map<std::string, Colour>, TextureMap> readMtl(const std::string &filename) {
	std::map<std::string, Colour> palette;
	TextureMap texture;
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
			int r = int(strtod(splitln[1].c_str(), nullptr) * 255);
			int g = int(strtod(splitln[2].c_str(), nullptr) * 255);
			int b = int(strtod(splitln[3].c_str(), nullptr) * 255);
			palette.insert({name, Colour(name, r, g, b)});
		}
		if (splitln[0] == "Ks") {
			palette.insert({name, Colour(name, 0, 0, 0)});
			palette[name].mirrored = true;
		}
		if (splitln[0] == "map_Kd") {
			// find texture filename and create TextureMap object
			size_t parent = filename.find_last_of('/');
			std::string texFile = filename.substr(0, parent + 1) + splitln[1];
			texture = TextureMap(texFile);
		}
	}

	return {palette, texture};
}

std::pair<std::vector<ModelTriangle>, TextureMap> readObj(const std::string &filename, float modelScale = 0.35) {
	std::vector<ModelTriangle> triangles;
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> texturePoints;
	std::map<std::string, Colour> palette;
	TextureMap texture;
	Colour trigColour;
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
			std::pair<std::map<std::string, Colour>, TextureMap> mtl = readMtl(matFile);
			palette = mtl.first;
			texture = mtl.second;
		}
		if (splitln[0] == "usemtl") {
			// look up colour from palette when reading 'usemtl'
			trigColour = palette[splitln[1]];
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
				tIndices.push_back(vt[1]);
			}
			// convert indices to ints and convert to zero-index (ignoring first char 'f')
			int iv0 = int(strtol(vIndices[0].c_str(), nullptr, 10) - 1);
			int iv1 = int(strtol(vIndices[1].c_str(), nullptr, 10) - 1);
			int iv2 = int(strtol(vIndices[2].c_str(), nullptr, 10) - 1);
			int it0 = int(strtol(tIndices[0].c_str(), nullptr, 10) - 1);
			int it1 = int(strtol(tIndices[1].c_str(), nullptr, 10) - 1);
			int it2 = int(strtol(tIndices[2].c_str(), nullptr, 10) - 1);
			// set vertices and their indices
			Vertex v0 = vertices[iv0], v1 = vertices[iv1], v2 = vertices[iv2];
			v0.index = iv0; v1.index = iv1; v2.index = iv2;
			// set texture points if used
			if(it0 != -1 && it1 != -1 && it2 != -1) {
				TexturePoint t0 = texturePoints[it0], t1 = texturePoints[it1], t2 = texturePoints[it2];
				v0.texturePoint = t0; v1.texturePoint = t1; v2.texturePoint = t2;
			}
			// calculate triangle face normal
			glm::vec3 normal = normalize(cross(v1.position - v0.position, v2.position - v0.position));
			// create new triangle
			ModelTriangle modelTriangle = {v0, v1, v2, trigColour};
			modelTriangle.normal = normal;
			modelTriangle.texture = texture.name;
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

	return {triangles, texture};
}


[[noreturn]] int main(int argc, char *argv[]) {
	// seed rand
	srand(time(nullptr));

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	// READ OBJ AND MTL FILES
	// create vector of triangles and textures
	std::vector<ModelTriangle> triangles;
	std::map<std::string, TextureMap> textures;

	// read obj files
	std::vector<ModelTriangle> box = readObj(BOX_OBJ).first;
	std::vector<ModelTriangle> sphere = readObj(SPH_OBJ).first;

	std::pair<std::vector<ModelTriangle>, TextureMap> hackspace = readObj(HACK_OBJ, 0.001f);
	std::vector<ModelTriangle> hackspaceLogo = hackspace.first;
	TextureMap hackTexture = hackspace.second;

	// insert obj triangles into vector of all triangles (if needed)
	triangles.insert(triangles.end(), box.begin(), box.end());
	triangles.insert(triangles.end(), sphere.begin(), sphere.end());
	triangles.insert(triangles.end(), hackspaceLogo.begin(), hackspaceLogo.end());

	// insert obj textures into vector of all textures
	textures.insert({hackTexture.name, hackTexture});

	// initialise light source and camera
	Camera camera = Camera(3.0f);
	PointLight light = PointLight();
	QuadLight ceiling = QuadLight();

	std::vector<std::reference_wrapper<Light>> lights;
	lights.emplace_back(ceiling);
	QuadLight quad = lights.at(0);
	std::cout << to_string(quad.e1) << std::endl; // doesn't get QuadLight values

	// setup fps counter
	time_t start = time(nullptr);
	time_t elapsed;
	int frames = 0;
	while (true) {
		// basic fps counter
		elapsed = time(nullptr) - start;
		if (elapsed >= 1) {
			std::cout << "fps: " << frames / elapsed << std::endl;
			start = time(nullptr);
			frames = 0;
		} frames++;
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, lights, camera, window);
		draw(triangles, textures, lights, camera, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}