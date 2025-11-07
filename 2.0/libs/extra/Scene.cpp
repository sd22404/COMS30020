#include "Scene.h"

static glm::vec3 vertexNormal(const Vertex &vertex, const std::vector<ModelTriangle> &triangles) {
	// for each vertex on triangles in the model, if it's the same as given vertex, add to list
	glm::vec3 total(0);
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

void Scene::moveLight(const Direction dir) const {
	switch (dir) {
		case UP:
			lights.at(0).position.y += 0.01;
			break;
		case DOWN:
			lights.at(0).position.y -= 0.01;
			break;
		case LEFT:
			lights.at(0).position.x -= 0.01;
			break;
		case RIGHT:
			lights.at(0).position.x += 0.01;
			break;
		case FORWARD:
			lights.at(0).position.z -= 0.01;
			break;
		case BACKWARD:
			lights.at(0).position.z += 0.01;
			break;
		default:
			break;
	}
}

Scene::Scene(const std::string &objFilename, std::vector<Light> &lights, const float modelScale) : lights(lights) {
	readObj(objFilename, modelScale);
	const auto it = textures.find("background");
	background = (it != textures.end()) ? &it->second : nullptr;
}

glm::vec3 Scene::backgroundColour(const float x, const float y) const {
	if (background == nullptr) return glm::vec3(0);
	return unpackColour(background->pixels.at(std::lround(y) * background->width + std::lround(x)));
}

RayTriangleIntersection Scene::closestIntersection(const Ray &ray, const float minDist, const float maxDist) const {
	float inverseClosestDistance = 0;
	int index = -1;
	RayTriangleIntersection intersection;
	intersection.triangleIndex = index;
	for (auto &triangle : triangles) {
		index++;
		// calculate edge vectors
		glm::vec3 e0 = triangle[1].position - triangle[0].position;
		glm::vec3 e1 = triangle[2].position - triangle[0].position;
		// calculate vector from startPoint to triangle
		glm::vec3 SPVector = ray.start - triangle[0].position;
		// generate direction/edge matrix
		glm::mat3 DEMatrix(-ray.dir, e0, e1);
		// find possible solution in [t, u, v]
		const glm::vec3 possibleSolution = inverse(DEMatrix) * SPVector;
		const float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
		// if closer than previously found solution, and within the bounds of the triangle, set new closest intersection
		if (t > minDist && t < maxDist && 1 / t > inverseClosestDistance && u >= 0 && u <= 1.0 && v >= 0 && v <= 1.0 && (u + v) <= 1.0) {
			intersection = RayTriangleIntersection(ray.start + t * ray.dir, t, triangle, index);
			intersection.u = u; intersection.v = v;
			inverseClosestDistance = 1 / t;
		}
	}

	auto &triangle = intersection.intersectedTriangle;
	const float u = intersection.u, v = intersection.v, w = 1.0f - u - v;

	if (!triangle.texture.empty()) {
		const TextureMap &texture = textures.find(triangle.texture)->second;
		const float texX = u * triangle[0].texturePoint.x + v * triangle[1].texturePoint.x + w * triangle[2].texturePoint.x;
		const float texY = u * triangle[0].texturePoint.y + v * triangle[1].texturePoint.y + w * triangle[2].texturePoint.y;
		const int iTexX = std::floor(texX * static_cast<float>(texture.width - 1));
		const int iTexY = std::floor(texY * static_cast<float>(texture.height - 1));
		intersection.textureSample = texture.pixels[iTexY * texture.width + iTexX];
	}

	if (!triangle.normalMap.empty()) {
		const TextureMap &texture = normalMaps.find(triangle.normalMap)->second;
		const float texX = u * triangle[0].normalPoint.x + v * triangle[1].normalPoint.x + w * triangle[2].normalPoint.x;
		const float texY = u * triangle[0].normalPoint.y + v * triangle[1].normalPoint.y + w * triangle[2].normalPoint.y;
		const int iTexX = std::floor(texX * static_cast<float>(texture.width - 1));
		const int iTexY = std::floor(texY * static_cast<float>(texture.height - 1));
		intersection.normalSample = texture.pixels[iTexY * texture.width + iTexX];
	}

	return intersection;
}


void Scene::readObj(const std::string &objFilename, float modelScale) {
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> texturePoints;
	TextureMap texture, normalMap;
	Material material;
	std::string line;
	std::vector<std::string> splitLn;
	std::ifstream file(objFilename);
    std::cout << "Reading model from " << objFilename << " ..." << std::endl;
	// read from file
	while(getline(file, line)) {
		// split line into segments
		splitLn = split(line, ' ');
		// read from mtl file
		if (splitLn[0] == "mtllib") {
			size_t parent = objFilename.find_last_of('/');
			std::string matFile = objFilename.substr(0, parent + 1) + splitLn[1];
			readMtl(matFile);
		}
		if (splitLn[0] == "usemtl") {
			// look up colour from palette
			material = materials[splitLn[1]];
			// check if texture and normal map exist for this material, otherwise reset them
			auto it = textures.find(splitLn[1]);
			if (it != textures.end()) texture = it->second;
			else texture = TextureMap();
			it = normalMaps.find(splitLn[1]);
			if (it != normalMaps.end()) normalMap = it->second;
			else normalMap = TextureMap();
		}
		if (splitLn[0] == "v") {
			// convert coords to floats and scale by 'modelScale' (ignoring first index which will be 'v')
			glm::vec3 vertex = modelScale * glm::vec3(strtod(splitLn[1].c_str(), nullptr), strtod(splitLn[2].c_str(), nullptr), strtod(splitLn[3].c_str(), nullptr));
			vertices.push_back(vertex);
		}
		if (splitLn[0] == "vt") {
			// convert to floats and create texture point (ignoring first index which will be 'vt')
			TexturePoint texturePoint = {std::stof(splitLn[1]), std::stof(splitLn[2])};
			texturePoints.push_back(texturePoint);
		}
		if (splitLn[0] == "f") {
			// split line into indices
			std::vector<std::string> vIndices, tIndices;
			for (size_t i = 1; i < splitLn.size(); i++) {
				std::vector<std::string> vt = split(splitLn[i], '/');
				if (!vt[0].empty()) vIndices.push_back(vt[0]);
				if (vt.size() > 1 && !vt[1].empty()) tIndices.push_back(vt[1]);
			}
			// convert indices to ints and zero-index (ignoring first char 'f')
			int iv0 = std::stoi(vIndices[0]) - 1;
			int iv1 = std::stoi(vIndices[1]) - 1;
			int iv2 = std::stoi(vIndices[2]) - 1;
			// set vertices and their indices
			auto v0 = Vertex(vertices[iv0]), v1 = Vertex(vertices[iv1]), v2 = Vertex(vertices[iv2]);
			v0.index = iv0; v1.index = iv1; v2.index = iv2;

			// get vertex texture indices
			if (tIndices.size() >= 3) {
				int it0 = std::stoi(tIndices[0]) - 1;
				int it1 = std::stoi(tIndices[1]) - 1;
				int it2 = std::stoi(tIndices[2]) - 1;
				TexturePoint t0 = texturePoints[it0], t1 = texturePoints[it1], t2 = texturePoints[it2];
				// set texture points if used
				if (!texture.name.empty()) {
					v0.texturePoint = t0; v1.texturePoint = t1; v2.texturePoint = t2;
				}
				// set normal map points if used
				if (!normalMap.name.empty()) {
					v0.normalPoint = t0; v1.normalPoint = t1; v2.normalPoint = t2;
				}
			}

			// calculate triangle face normal
			glm::vec3 normal = normalize(cross(v1.position - v0.position, v2.position - v0.position));
			// create new triangle
			auto modelTriangle = ModelTriangle(v0, v1, v2, material);
			modelTriangle.normal = normal;
			modelTriangle.texture = texture.name.empty() ? "" : material.name;
			modelTriangle.normalMap = normalMap.name.empty() ? "" : material.name;
			// add triangle
			triangles.emplace_back(modelTriangle);
		}
	}

    file.close();

	// calculate vertex normals based on neighbouring triangles
	for (auto &triangle : triangles) {
		for (auto &v : triangle.vertices) {
			v.normal = vertexNormal(v, triangles);
		}
	}

    std::cout << "Read " << triangles.size() << " triangles from " << objFilename << std::endl;
}

void Scene::readMtl(const std::string &filename) {
	std::ifstream file(filename);
	std::string line;
	std::string name;
	std::vector<std::string> splitLn;
	// read from file
	while(getline(file, line)) {
		splitLn = split(line, ' ');
		if (splitLn[0] == "newmtl") {
			// set colour name
			name = splitLn[1];
		}
		if (splitLn[0] == "Kd") {
			// create new colour from given values
			float r = strtof(splitLn[1].c_str(), nullptr);
			float g = strtof(splitLn[2].c_str(), nullptr);
			float b = strtof(splitLn[3].c_str(), nullptr);
			materials.insert({name, Material(name, glm::vec3(r, g, b))});
		}
		if (splitLn[0] == "Ks") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			materials[name].reflectivity = strtof(splitLn[1].c_str(), nullptr);
		}
		if (splitLn[0] == "Ke") {
			materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
			materials[name].emissive = true;
		}
		if (splitLn[0] == "Ka") {
			materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
			materials[name].glassy = true;
			materials[name].refractiveIndex = strtof(splitLn[1].c_str(), nullptr);
		}
		if (splitLn[0] == "map_Kd") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			// find texture filename and create TextureMap object
			/*
			size_t parent = filename.find_last_of('/');
			std::string texFile = filename.substr(0, parent + 1) + splitLn[1];
			*/
			std::string texFile = "../assets/textures/" + splitLn[1];
			auto texture = TextureMap(texFile, name);
			textures.insert({name, texture});
		}
		if (splitLn[0] == "map_bump") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			// find texture filename and create TextureMap object
			/*size_t parent = filename.find_last_of('/');
			std::string texFile = filename.substr(0, parent + 1) + splitLn[1];
			*/
			std::string texFile = "../assets/normals/" + splitLn[1];
			auto normal = TextureMap(texFile, name);
			normalMaps.insert({name, normal});
		}
	}

    file.close();

    std::cout << "Read " << materials.size() << " materials from " << filename << std::endl;
}
