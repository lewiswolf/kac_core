// core
#include <stdlib.h>	 // srand
#include <time.h>	 // time

#include <iostream>	 // cout, endl

// src
#include "../src/geometry.hpp"

int main() {
	// init random number generator
	srand(time(NULL));

	std::cout << "Test generateConvexPolygon: ";
	Vertices v = generateConvexPolygon(7);
	for (int i = 0; i < v.size(); i++) {
		std::cout << "(" << v[i].x << ", " << v[i].y << ") ";
	}
	std::cout << std::endl;

	std::cout << "isConvex: ";
	std::cout << (isConvex(v) ? "true" : "false");
	std::cout << std::endl << std::endl;

	return 0;
}