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
	const int efficiency = 1000;
	std::cout << "\nEfficiency for " << efficiency << " vertices...\n";
	Vertices test_v;
	{
		Timer timer("	generateConvexPolygon");
		test_v = generateConvexPolygon(efficiency);
	}
	{
		Timer timer("	isConvex");
		isConvex(test_v);
	}
	{
		Timer timer("	isColinear");
		isColinear(test_v[efficiency - 1], test_v[0], test_v[1]);
		for (unsigned int i = 1; i < efficiency - 1; i++) {
			isColinear(test_v[i - 1], test_v[i], test_v[i + 1]);
		};
		isColinear(test_v[efficiency - 2], test_v[efficiency - 1], test_v[0]);
	}
}