#include "RenderUtils.h"

std::vector<CanvasPoint> RenderUtils::interpolateCanvasPoints(CanvasPoint from, CanvasPoint to, float numberOfValues) {
	// if one or fewer values, return only the start point (as it will be the same as the end point)
	std::vector<CanvasPoint> result = {from};
	if (numberOfValues < 1) return result;
	CanvasPoint interval = {to.x - from.x, to.y - from.y, to.depth - from.depth};
	float stepX = interval.x / numberOfValues;
	float stepY = interval.y / numberOfValues;
	float stepD = interval.depth / numberOfValues;

	std::vector<TexturePoint> tps = interpolateTexturePoints(from.texturePoint, to.texturePoint, numberOfValues);
	for (float i = 1; i < numberOfValues; i++) {
		CanvasPoint point = {from.x + i * stepX, from.y + i * stepY, from.depth + i * stepD};
		point.texturePoint = tps[i];
		result.push_back(point);
	}
	return result;
}

std::vector<TexturePoint> RenderUtils::interpolateTexturePoints(TexturePoint from, TexturePoint to, float numberOfValues) {
	// if one or fewer values, return only the start point (as it will be the same as the end point)
	std::vector<TexturePoint> result = {from};
	if (numberOfValues < 1) return result;

	TexturePoint interval = {to.x - from.x, to.y - from.y};
	float stepX = interval.x / numberOfValues;
	float stepY = interval.y / numberOfValues;

	for (float i = 1; i < numberOfValues; i++) {
		TexturePoint point = {from.x + i * stepX, from.y + i * stepY};
		result.push_back(point);
	}
	return result;
}