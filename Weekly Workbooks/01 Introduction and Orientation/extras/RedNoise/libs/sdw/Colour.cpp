#include "Colour.h"
#include <utility>

Colour::Colour() = default;
Colour::Colour(int r, int g, int b) : red(r), green(g), blue(b) {}
Colour::Colour(std::string n, int r, int g, int b) :
		name(std::move(n)),
		red(r), green(g), blue(b) {}

Colour operator*(const float scale, const Colour &col) {
	Colour newCol = Colour(col);
	newCol.red *= scale;
	newCol.green *= scale;
	newCol.blue *= scale;
	return newCol;
}

Colour operator*(const Colour &col1, const Colour &col2) {
	Colour newCol = Colour(col1);
	newCol.red *= col2.red;
	newCol.green *= col2.green;
	newCol.blue *= col2.blue;
	return newCol;
}

std::ostream &operator<<(std::ostream &os, const Colour &colour) {
	os << colour.name << " ["
	   << colour.red << ", "
	   << colour.green << ", "
	   << colour.blue << "]";
	return os;
}