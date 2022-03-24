// core
#include <iostream>
#include <stdlib.h>

// src
#include "../../src/geometry.hpp"

// test
#include "./utils.hpp"

void testShapes() {
	// tests
	int N = 10;
	Vertices v = generateConvexPolygon(N);
	std::cout << "Test if generatedConvexPolygon produces n vertices... ";
	booleanTest(v.size() == N);
	std::cout << "Test if generatedConvexPolygon is convex... ";
	booleanTest(isConvex(v));

	// efficiency
	std::cout << "\nEfficiency for 1000 vertices...\n";
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
}