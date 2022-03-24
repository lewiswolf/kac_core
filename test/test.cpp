// core
#include <iostream>	 // cout, endl
#include <stdlib.h>	 // srand
#include <time.h>	 // time

// src
#include "../src/geometry.hpp"

// test
#include "./src/timer.hpp"

int main() {
	// init random number generator
	srand(time(NULL));

	std::cout << "Test generateConvexPolygon: ";
	Vertices v = generateConvexPolygon(7);

	for (int i = 0; i < v.size(); i++) {
		std::cout << "(" << v[i].x << ", " << v[i].y << ") ";
	}
	std::cout << "isConvex: ";
	std::cout << (isConvex(v) ? "true" : "false") << "\n\n";

	std::cout << "for 1000 vertices...\n";
	Vertices test;
	{
		Timer timer("	generateConvexPolygon");
		test = generateConvexPolygon(1000);
	}
	{
		Timer timer("	isConvex");
		isConvex(test);
	}
	{
		Timer timer("	isColinear");
		isColinear(test[999], test[0], test[1]);
		for (unsigned int i = 1; i < 999; i++) {
			isColinear(test[i - 1], test[i], test[i + 1]);
		};
		isColinear(test[998], test[999], test[0]);
	}

	// std::cout << "bessel J_0(4): ";
	// std::cout << bessel(0, 4);
	// std::cout << std::endl << std::endl;

	return 0;
}