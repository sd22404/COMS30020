#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Utils.h"
#include <cstdint>

struct TextureMap {
	std::string name;
	size_t width;
	size_t height;
	std::vector<uint32_t> pixels;

	TextureMap();
	TextureMap(const std::string &filename, const std::string &name);
	friend std::ostream &operator<<(std::ostream &os, const TextureMap &point);
};
