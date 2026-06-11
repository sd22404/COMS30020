#include "HitInfo.h"
#include <utility>

HitInfo::HitInfo() = default;

HitInfo::HitInfo(const glm::vec3 &point, float distance, int modelIndex, int triangleIndex) :
	point(point),
	distance(distance),
	u(0), v(0), w(0), modelIndex(modelIndex), triIndex(triangleIndex) {}

std::ostream &operator<<(std::ostream &os, const HitInfo &intersection) {
	os << "Intersection is at [" << intersection.point[0] << "," << intersection.point[1] << "," <<
	   intersection.point[2] << "] on model " << intersection.modelIndex <<
	   " at a distance of " << intersection.distance;
	return os;
}
