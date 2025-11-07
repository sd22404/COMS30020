#include "RayTriangleIntersection.h"

#include <utility>

RayTriangleIntersection::RayTriangleIntersection() = default;
RayTriangleIntersection::RayTriangleIntersection(const glm::vec3 &point, float distance, ModelTriangle triangle, size_t index) :
	intersectionPoint(point),
	distanceFromCamera(distance),
	intersectedTriangle(std::move(triangle)),
	triangleIndex(index), u(0), v(0), textureSample(0), normalSample(0) {}

std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection) {
	os << "Intersection is at [" << intersection.intersectionPoint[0] << "," << intersection.intersectionPoint[1] << "," <<
	   intersection.intersectionPoint[2] << "] on triangle " << intersection.triangleIndex <<
	   " at a distance of " << intersection.distanceFromCamera;
	return os;
}
