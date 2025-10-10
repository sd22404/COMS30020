#include "RenderUtils.h"

std::vector<CanvasPoint> RenderUtils::interpolateCanvasPoints(CanvasPoint from, CanvasPoint to, size_t numberOfValues) {
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

std::vector<TexturePoint> RenderUtils::interpolateTexturePoints(TexturePoint from, TexturePoint to, size_t numberOfValues) {
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