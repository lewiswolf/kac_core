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
	// geometry/generate_polygon.hpp
	std::cout << "\nProfiler for `./geometry/generate_polygon.hpp`.\n";
	std::cout << "Efficiency relative to " << N << " vertices...\n";
	T::Polygon P_convex;
	T::Polygon P;
	T::Polygon P_tmp;
	{
		Timer timer("  generateConvexPolygon");
		P_convex = g::generateConvexPolygon(N);
	}
	{
		Timer timer("  generateIrregularStar");
		P_tmp = g::generateIrregularStar(N);
	}
	{
		Timer timer("  generatePolygon");
		P = g::generatePolygon(N);
	}
	{
		Timer timer("  generateUnitRectangle");
		P_tmp = g::generateUnitRectangle(0.5);
	}

	// geometry/morphisms.hpp
	std::cout << "\nProfiler for `./geometry/morphisms.hpp`.\n";
	std::cout << "Efficiency relative to a " << N << " sided polygon...\n";
	{
		Timer timer("  normalisePolygon");
		g::normalisePolygon(P);
	}
	{
		Timer timer("  normaliseConvexPolygon");
		g::normaliseConvexPolygon(P_convex);
	}
	{
		Timer timer("  normaliseSimplePolygon");
		g::normaliseSimplePolygon(P);
	}

	// geometry/polygon_properties.hpp
	std::cout << "\nProfiler for `./geometry/polygon_properties.hpp`.\n";
	std::cout << "Efficiency relative to a " << N << " sided polygon...\n";
	T::Point centroid = g::polygonCentroid(P);
	{
		Timer timer("  isConvex");
		g::isConvex(P_convex);
	}
	{
		Timer timer("  isPointInsideConvexPolygon");
		g::isPointInsideConvexPolygon(centroid, P_convex);
	}
	{
		Timer timer("  isSimple");
		g::isSimple(P);
	}
	{
		Timer timer("  largestVector");
		g::largestVector(P);
	}
	{
		Timer timer("  polygonCentroid");
		g::polygonCentroid(P);
	}
	{
		Timer timer("  polygonArea");
		g::polygonArea(P);
	}

	return 0;
}
