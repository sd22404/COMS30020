#include "Scene.h"

static glm::vec3 vertexNormal(Vertex &vertex, std::vector<ModelTriangle> &triangles) {
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

void Scene::moveLight(Direction dir) {
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

Scene::Scene(const std::string &objFilename, std::vector<Light> &lights, float modelScale) : lights(lights) {
	readObj(objFilename, modelScale);
	auto it = textures.find("background");
	background = (it != textures.end()) ? &it->second : nullptr;
}

uint32_t Scene::backgroundColour(float x, float y) {
	if (background == nullptr) return 0;
	return background->pixels.at((y * background->width) + x);
}

RayTriangleIntersection Scene::closestIntersection(Ray &ray, float minDist, float maxDist) {
	float inverseClosestDistance = 0;
	int index = -1;
	RayTriangleIntersection intersection;
	intersection.triangleIndex = index;
	for (auto &triangle : triangles) {
		index++;
		// calculate edge vectors
		glm::vec3 e0 = triangle.v1().position - triangle.v0().position;
		glm::vec3 e1 = triangle.v2().position - triangle.v0().position;
		// calculate vector from startPoint to triangle
		glm::vec3 SPVector = ray.start - triangle.v0().position;
		// generate direction/edge matrix
		glm::mat3 DEMatrix(-ray.dir, e0, e1);
		// find possible solution in [t, u, v]
		glm::vec3 possibleSolution = inverse(DEMatrix) * SPVector;
		float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
		// if closer than previously found solution, and within the bounds of the triangle, set new closest intersection
		if (t > minDist && t < maxDist && 1 / t > inverseClosestDistance && u >= 0 && u <= 1.0 && v >= 0 && v <= 1.0 && (u + v) <= 1.0) {
			intersection = RayTriangleIntersection(ray.start + t * ray.dir, t, triangle, index);
			intersection.u = u; intersection.v = v;
			inverseClosestDistance = 1 / t;
		}
	}

	return intersection;
}


void Scene::readObj(const std::string &filename, float modelScale) {
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> texturePoints;
	TextureMap texture, normalMap;
	Material material;
	std::string line;
	std::vector<std::string> splitln;
	std::ifstream file(filename);
    std::cout << "Reading model from " << filename << " ..." << std::endl;
	// read from file
	while(getline(file, line)) {
		// split line into segments
		splitln = split(line, ' ');
		// read from mtl file
		if (splitln[0] == "mtllib") {
			size_t parent = filename.find_last_of('/');
			std::string matFile = filename.substr(0, parent + 1) + splitln[1];
			readMtl(matFile);
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
				if (!vt[0].empty()) vIndices.push_back(vt[0]);
				if (!vt[1].empty()) tIndices.push_back(vt[1]);
			}
			// convert indices to ints and zero-index (ignoring first char 'f')
			int iv0 = int(strtol(vIndices[0].c_str(), nullptr, 10) - 1);
			int iv1 = int(strtol(vIndices[1].c_str(), nullptr, 10) - 1);
			int iv2 = int(strtol(vIndices[2].c_str(), nullptr, 10) - 1);
			// set vertices and their indices
			Vertex v0 = Vertex(vertices[iv0]), v1 = Vertex(vertices[iv1]), v2 = Vertex(vertices[iv2]);
			v0.index = iv0; v1.index = iv1; v2.index = iv2;
			// set texture points if used
			if (tIndices.size() >= 3) {
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
			ModelTriangle modelTriangle = ModelTriangle(v0, v1, v2, material);
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

    std::cout << "Read " << triangles.size() << " triangles from " << filename << std::endl;
}

void Scene::readMtl(const std::string &filename) {
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
			normalMaps.insert({name, normal});
		}
	}

    file.close();

    std::cout << "Read " << materials.size() << " materials from " << filename << std::endl;
}
