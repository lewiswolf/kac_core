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
namespace p = kac_core::physics;	 // physics

// test
#include "./utils.hpp"

const unsigned long N = 200;
const unsigned long t = 48000;

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

	// geometry/mappings.hpp
	std::cout << "\nProfiler for `./geometry/mappings.hpp`.\n";
	std::cout << "Efficiency relative to " << N << " points...\n";
	{
		Timer timer("  circle2Square");
		for (unsigned long n = 0; n < P.size(); n++) { g::simpleElliptic_Circle2Square(P[n]); }
	}
	{
		Timer timer("  square2Circle");
		for (unsigned long n = 0; n < P.size(); n++) { g::simpleElliptic_Square2Circle(P[n]); }
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
	{
		Timer timer("  scalePolygonByArea");
		g::scalePolygonByArea(P, 100.0);
	}

	// geometry/polygon_properties.hpp
	std::cout << "\nProfiler for `./geometry/polygon_properties.hpp`.\n";
	std::cout << "Efficiency relative to a " << N << " sided polygon...\n";
	T::Point centroid = g::polygonCentroid(P);
	T::Point convex_centroid = g::polygonCentroid(P_convex);
	{
		Timer timer("  isConvex");
		g::isConvex(P_convex);
	}
	{
		Timer timer("  isPointInsideConvexPolygon");
		g::isPointInsideConvexPolygon(convex_centroid, P_convex);
	}
	{
		Timer timer("  isPointInsidePolygon");
		g::isPointInsidePolygon(centroid, P);
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

	// ./physics/modes
	std::cout << "\nProfiler for `./physics/modes.\n";
	std::cout << "Efficiency relative to " << N << " modes and a waveform of duration of " << t
			  << " samples...\n";
	T::Matrix_1D F_1d = p::linearSeries(N);
	T::Matrix_1D A_1d = p::linearAmplitudes(0.5, N);
	{
		Timer timer("  WaveEquationWaveform2D");
		T::Matrix_1D waveform = p::WaveEquationWaveform1D(F_1d, A_1d, 1.0, 1 / t, t);
	}
	std::cout << "Efficiency relative to a " << N << " X " << N << " matrix of modes...\n";
	{
		Timer timer("  circularChladniPattern");
		T::BooleanImage circular_pattern = p::circularChladniPattern(2, 2, N, 0.1);
	}
	{
		Timer timer("  circularCymatics");
		T::Matrix_2D circular_pattern = p::circularCymatics(2, 2, N);
	}
	{
		Timer timer("  rectangularChladniPattern");
		T::BooleanImage rectangular_pattern = p::rectangularChladniPattern(2, 2, N, N, 0.1);
	}
	{
		Timer timer("  rectangularCymatics");
		T::Matrix_2D rectangular_pattern = p::rectangularCymatics(2, 2, N, N);
	}
	std::cout << "Efficiency relative to a " << N << " X " << N
			  << " matrix of modes and a waveform of duration of " << t << " samples...\n";
	T::Matrix_2D F_2d = p::rectangularSeries(N, N, 1);
	T::Matrix_2D A_2d = p::rectangularAmplitudes(0.5, 0.5, N, N, 1);
	{
		Timer timer("  WaveEquationWaveform2D");
		T::Matrix_1D waveform = p::WaveEquationWaveform2D(F_2d, A_2d, 1.0, 1 / t, t);
	}

	return 0;
}
