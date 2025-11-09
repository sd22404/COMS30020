#include "HitInfo.h"
#include <utility>

HitInfo::HitInfo() = default;

HitInfo::HitInfo(const glm::vec3 &point, float distance, int modelIndex, int triangleIndex) :
	intersectionPoint(point),
	distanceFromCamera(distance),
	u(0), v(0), w(0), modelIndex(modelIndex), triIndex(triangleIndex) {}

std::ostream &operator<<(std::ostream &os, const HitInfo &intersection) {
	os << "Intersection is at [" << intersection.intersectionPoint[0] << "," << intersection.intersectionPoint[1] << "," <<
	   intersection.intersectionPoint[2] << "] on model " << intersection.modelIndex <<
	   " at a distance of " << intersection.distanceFromCamera;
	return os;
}
