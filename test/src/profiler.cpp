/*
Tests and profiling for /shapes.
*/

// core
#include <iostream>
#include <numbers>
#include <stdlib.h>
#include <string>

// src
#include <kac_core.hpp>
namespace T = kac_core::types;		 // types
namespace g = kac_core::geometry;	 // geometry
namespace p = kac_core::physics;	 // physics

// test
#include "./utils.hpp"

// magic numbers
const std::size_t N = 200;
const std::size_t t = 48000;
// implicit magic
const std::string N_string = " " + std::to_string(N) + " ";
const std::string t_string = " " + std::to_string(t) + " ";

int main() {
	// geometry/generate_polygon.hpp
	printColouredText("\nProfiler for `./geometry/generate_polygon.hpp`.", 36);
	printColouredText("Efficiency relative to" + N_string + "vertices...", 35);
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
		Timer timer("  generateRegularPolygon");
		P_tmp = g::generateRegularPolygon(N);
	}
	{
		Timer timer("  generateUnitRectangle");
		P_tmp = g::generateUnitRectangle(0.5);
	}
	{
		Timer timer("  generateUnitTriangle");
		P_tmp = g::generateUnitTriangle(1., std::numbers::pi * 0.5);
	}

	// geometry/mappings.hpp
	printColouredText("\nProfiler for `./geometry/mappings.hpp`.", 36);
	printColouredText("Efficiency relative to" + N_string + "points...", 35);
	{
		Timer timer("  circleToSquare");
		for (std::size_t n = 0; n < P.size(); n++) { g::circleToSquare(P[n]); }
	}
	{
		Timer timer("  squareToCircle");
		for (std::size_t n = 0; n < P.size(); n++) { g::squareToCircle(P[n]); }
	}
	{
		Timer timer("  squareToTriangle");
		for (std::size_t n = 0; n < P.size(); n++) { g::squareToTriangle(P[n]); }
	}
	{
		Timer timer("  triangleToSquare");
		for (std::size_t n = 0; n < P.size(); n++) { g::triangleToSquare(P[n]); }
	}

	// geometry/morphisms.hpp
	printColouredText("\nProfiler for `./geometry/morphisms.hpp`.", 36);
	printColouredText("Efficiency relative to a" + N_string + "sided polygon...", 35);
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
		g::scalePolygonByArea(P, 100.);
	}

	// geometry/polygon_properties.hpp
	printColouredText("\nProfiler for `./geometry/polygon_properties.hpp`.", 36);
	printColouredText("Efficiency relative to a" + N_string + "sided polygon...", 35);
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
	printColouredText("\nProfiler for `./physics/modes`.", 36);
	printColouredText("Efficiency relative to" + N_string + "X" + N_string + "modes...", 35);
	{
		Timer timer("  linearAmplitudes");
		T::Matrix_1D linear_pattern = p::linearAmplitudes(0.5, N);
	}
	{
		Timer timer("  linearCymatics");
		T::Matrix_1D linear_pattern = p::linearCymatics(2, N);
	}
	{
		Timer timer("  linearSeries");
		T::Matrix_1D linear_pattern = p::linearSeries(N);
	}
	T::Matrix_2D _S = p::circularSeries(N, N);
	{
		Timer timer("  circularAmplitudes");
		T::Matrix_2D circular_pattern = p::circularAmplitudes(0.5, 0.5, _S);
	}
	{
		Timer timer("  circularCymatics");
		T::Matrix_2D circular_pattern = p::circularCymatics(2, 2, N);
	}
	{
		Timer timer("  circularSeries");
		T::Matrix_2D circular_pattern = p::circularSeries(N, N);
	}
	{
		Timer timer("  rectangularAmplitudes");
		T::Matrix_2D rectangular_pattern = p::rectangularAmplitudes(0.5, 0.5, N, N);
	}
	{
		Timer timer("  rectangularCymatics");
		T::Matrix_2D rectangular_pattern = p::rectangularCymatics(2., 2., N, N);
	}
	{
		Timer timer("  rectangularSeries");
		T::Matrix_2D rectangular_pattern = p::rectangularSeries(N, N, 1.);
	}
	printColouredText(
		"Efficiency relative to" + N_string + "modes and a waveform" + t_string
			+ "samples in length...",
		35
	);
	T::Matrix_1D F_1d = p::linearSeries(N);
	T::Matrix_1D A_1d = p::linearAmplitudes(0.5, N);
	{
		Timer timer("  AdditiveSynthesis1D");
		T::Matrix_1D waveform = p::AdditiveSynthesis1D(F_1d, A_1d, 1., 1 / t, t);
	}
	printColouredText(
		"Efficiency relative to" + N_string + "X" + N_string + "modes and a waveform" + t_string
			+ "samples in length...",
		35
	);
	T::Matrix_2D F_2d = p::rectangularSeries(N, N, 1.);
	T::Matrix_2D A_2d = p::rectangularAmplitudes(0.5, 0.5, N, N);
	{
		Timer timer("  AdditiveSynthesis2D");
		T::Matrix_1D waveform = p::AdditiveSynthesis2D(F_2d, A_2d, 1., 1 / t, t);
	}

	// ./physics/fdtd
	printColouredText("\nProfiler for `./physics/fdtd`.", 36);
	printColouredText("Efficiency relative to an" + N_string + "X" + N_string + "matrix...", 35);
	{
		Timer timer("  raisedCosine1D");
		T::Matrix_1D waveform = p::raisedCosine1D(0.5, 0.1, N);
	}
	{
		Timer timer("  raisedCosine2D");
		T::Matrix_2D waveform = p::raisedCosine2D(T::Point(0.5, 0.5), 0.1, N, N);
	}
	{
		Timer timer("  raisedTriangle1D");
		T::Matrix_1D waveform = p::raisedTriangle1D(0.5, 0.1, 0.1, N);
	}
	{
		Timer timer("  raisedTriangle2D");
		T::Matrix_2D waveform = p::raisedTriangle2D(T::Point(0.5, 0.5), 0.1, 0.1, 0.1, 0.1, N, N);
	}
	printColouredText(
		"Efficiency relative to a" + N_string + "matrix simulation and a waveform" + t_string
			+ "samples in length...",
		35
	);
	T::Matrix_1D u1_1(N, 0.);
	T::Matrix_1D u1_0(N, 0.);
	u1_1[(std::size_t)N * 0.5] = 1.;
	{
		Timer timer("  FDTDWaveform1D");
		p::FDTDWaveform1D(u1_0, u1_1, 1., 2. - 4., 1., t, 0.5);
	}
	printColouredText(
		"Efficiency relative to a" + N_string + "X" + N_string + "matrix simulation and a waveform"
			+ t_string + "samples in length...",
		35
	);
	double cfl_2 = pow(1 / std::numbers::sqrt2, 2.);
	T::Matrix_2D u2_1(N, std::vector<double>(N, 0.));
	T::Matrix_2D u2_0(N, std::vector<double>(N, 0.));
	u2_1[(std::size_t)N * 0.5][(std::size_t)N * 0.5] = 1.;
	T::BooleanImage_2D B(N, std::vector<short>(N, 1));
	// dirichlet boundary
	for (std::size_t i = 0; i < N; i++) { B[i][0] = B[i][N - 1] = 0; }
	for (std::size_t j = 0; j < N; j++) { B[0][j] = B[N - 1][j] = 0; }
	{
		Timer timer("  FDTDWaveform2D");
		p::FDTDWaveform2D(u2_0, u2_1, B, cfl_2, 2. - 4. * cfl_2, 1., t, T::Point(0.5, 0.5));
	}

	return 0;
}
