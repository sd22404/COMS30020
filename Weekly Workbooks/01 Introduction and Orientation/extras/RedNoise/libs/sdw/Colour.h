#pragma once

#include <iostream>

struct Colour {
	std::string name;
	int red{};
	int green{};
	int blue{};
	bool mirrored{};
	bool glassy{};
	float roughness{};
	Colour();
	Colour(int r, int g, int b);
	Colour(std::string n, int r, int g, int b);

	friend Colour operator*(const float scale, const Colour &col);
	friend Colour operator*(const Colour &col1, const Colour &col2);
};

std::ostream &operator<<(std::ostream &os, const Colour &colour);
