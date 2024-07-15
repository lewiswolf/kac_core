/*
Tests and profiling for /shapes.
*/

// core
#include <iostream>
#include <stdlib.h>

// src
#include <kac_core.hpp>
namespace T = kac_core::types;		 // types
namespace g = kac_core::geometry;	 // geometry

// test
#include "./utils.hpp"

unsigned long N = 200;

int main() {
	std::cout << "Profiler for `./src/shapes`.\n";
	std::cout << "Efficiency relative to " << N << " vertices...\n";
	T::Polygon P;

	// geometry/generate_polygon.hpp
	{
		Timer timer("  generateConvexPolygon");
		P = g::generateConvexPolygon(N);
	}
	{
		Timer timer("  generateIrregularStar");
		P = g::generateIrregularStar(N);
	}
	{
		Timer timer("  generatePolygon");
		P = g::generatePolygon(N);
	}

	return 0;
}
