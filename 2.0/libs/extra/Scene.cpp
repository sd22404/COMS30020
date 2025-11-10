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

Scene::Scene(const std::vector<Obj> &objs, std::vector<Light> &lights) : lights(lights) {
	for (const auto &obj : objs) {
		readObj(obj);
	}
}

glm::vec3 Scene::backgroundColour(const int x, const int y) const {
	if (background == nullptr) return glm::vec3(0);
	return unpackColour(background->pixels.at(y * background->width + x));
}

HitInfo Scene::closestIntersection(const Ray &ray, const float minDist, const float maxDist) const {
	float inverseClosestDistance = 0;
	int modelIndex = -1;
	HitInfo intersection;
	for (auto &model : models) {
		modelIndex++;
		int triIndex = -1;
		for (auto &triangle : model.triangles) {
			triIndex++;
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
			if (t > minDist && t < maxDist && 1.0f / t > inverseClosestDistance && u >= 0 && u <= 1.0 && v >= 0 && v <= 1.0 && (u + v) <= 1.0) {
				intersection = HitInfo(ray.start + t * ray.dir, t, modelIndex, triIndex);
				intersection.u = u; intersection.v = v; intersection.w = 1.0f - u - v;
				inverseClosestDistance = 1 / t;
			}
		}
	}

	return intersection;
}


void Scene::readObj(const Obj &obj) {
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> texturePoints;
	std::unordered_map<std::string, Material> materials;
	int count = 0;
	std::string line;
	std::vector<std::string> splitLn;
	std::ifstream file(obj.filename);
	// read from file
	while(getline(file, line)) {
		// split line into segments
		splitLn = split(line, ' ');
		// read from mtl file
		if (splitLn[0] == "mtllib") {
			size_t parent = obj.filename.find_last_of('/');
			std::string matFile = obj.filename.substr(0, parent + 1) + splitLn[1];
			materials = readMtl(matFile);
		}
		if (splitLn[0] == "o") {
			// Start a new model entry
			models.emplace_back();
			models.back().name = splitLn[1];
			models.back().sMode = obj.shadingMode;
			count++;
		}
		if (splitLn[0] == "usemtl") {
			// assign material to current model
			auto &dst = models.back().material;
			const auto &src = materials[splitLn[1]];
			dst.name = src.name;
			dst.diffuse = src.diffuse;
			dst.specular = src.specular;
			dst.shininess = src.shininess;
			dst.ambient = src.ambient;
			dst.reflectivity = src.reflectivity;
			dst.texture = src.texture;
			dst.normalMap = src.normalMap;
		}
		if (splitLn[0] == "v") {
			// convert coords to floats and scale by 'modelScale' (ignoring first index which will be 'v')
			glm::vec3 vertex = obj.offset + obj.scale * glm::vec3(stof(splitLn[1]), stof(splitLn[2]), stof(splitLn[3]));
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
				// set texture points
				v0.texturePoint = t0; v1.texturePoint = t1; v2.texturePoint = t2;
			}

			// calculate triangle face normal
			glm::vec3 normal = normalize(cross(v1.position - v0.position, v2.position - v0.position));
			// create new triangle
			auto modelTriangle = ModelTriangle(v0, v1, v2);
			modelTriangle.normal = normal;
			// add triangle
			models.back().triangles.emplace_back(modelTriangle);
		}
	}

    file.close();

	// calculate vertex normals based on neighbouring triangles
	for (auto &model : models) {
		for (auto &triangle : model.triangles) {
			for (auto &v : triangle.vertices) {
				v.normal = vertexNormal(v, model.triangles);
			}
		}
	}

    std::cout << "Read " << count << " model(s) from " << obj.filename << std::endl;
}

std::unordered_map<std::string, Material> Scene::readMtl(const std::string &filename) {
	std::ifstream file(filename);
	std::string line;
	std::vector<std::string> splitLn;
	std::string name;
	std::unordered_map<std::string, Material> materials;
	// read from file
	while(getline(file, line)) {
		splitLn = split(line, ' ');
		if (splitLn[0] == "newmtl") {
			// set colour name
			name = splitLn[1];
		}
		if (splitLn[0] == "Kd") {
			// create new colour from given values
			float r = stof(splitLn[1]);
			float g = stof(splitLn[2]);
			float b = stof(splitLn[3]);
			materials.insert({name, Material(name, glm::vec3(r, g, b))});
		}
		if (splitLn[0] == "Ks") {
			float r = stof(splitLn[1]);
			float g = stof(splitLn[2]);
			float b = stof(splitLn[3]);
			materials.insert({name, Material(name, glm::vec3(r, g, b))});
			materials[name].specular = glm::vec3(r, g, b);
		}
		if (splitLn[0] == "Ns") {
			materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
			materials[name].shininess = stof(splitLn[1]);
		}
		if (splitLn[0] == "d") {
			materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
			materials[name].reflectivity = stof(splitLn[1]);
		}
		if (splitLn[0] == "Ni") {
			materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
			materials[name].glassy = true;
			materials[name].refractiveIndex = stof(splitLn[1]);
		}
		// if (splitLn[0] == "Ke") {
		// 	materials.insert({name, Material(name, glm::vec3(1, 1, 1))});
		// 	materials[name].emissive = true;
		// }
		if (splitLn[0] == "map_Kd") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			std::string texFile = "./assets/textures/" + splitLn[1];
			materials[name].texture = TextureMap(texFile, name);
		}
		if (splitLn[0] == "map_bump") {
			materials.insert({name, Material(name, glm::vec3(0, 0, 0))});
			std::string texFile = "./assets/normals/" + splitLn[1];
			materials[name].normalMap = TextureMap(texFile, name);
		}
	}

    file.close();

    std::cout << "Read " << materials.size() << " material(s) from " << filename << std::endl;
	return materials;
}
