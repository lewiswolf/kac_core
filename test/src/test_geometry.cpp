/*
Tests for /geometry.
*/

#include <iostream>

// src
#include <kac_core.hpp>
namespace T = kac_core::types;		 // types
namespace g = kac_core::geometry;	 // geometry

// test
#include "./utils.hpp"

int main() {
	/*
	Initialise polygons.
	*/
	unsigned long N = 10;
	T::Polygon P_convex = g::generateConvexPolygon(N, 1);
	/*
	Test the seed argument of generateConvexPolygon.
	*/
	// T::Matrix_2D expected = {
	// 	{-0.552113, 0.713991},
	// 	{-0.670506, 0.425103},
	// 	{-0.713519, 0.302169},
	// 	{-0.777715, -0.0903147},
	// 	{-0.745129, -0.162445},
	// 	{-0.621568, -0.323642},
	// 	{0.0841093, -0.713991},
	// 	{0.777715, 0.331095},
	// 	{0.481532, 0.648682},
	// 	{0.335852, 0.677941}
	// };
	// batchBooleanTest(
	// 	"generatedConvexPolygon produces the expected output for a given seed",
	// 	10,
	// 	[&P_convex, &expected](unsigned int n) {
	// 		return abs(P_convex[n].x - expected[n][0]) < 0.01
	// 			&& abs(P_convex[n].y - expected[n][1]) < 0.01;
	// 	}
	// );
	/*
	Test the properties of generateConvexPolygon.
	*/
	booleanTest("generatedConvexPolygon produces n vertices", P_convex.size() == N);
	booleanTest("generatedConvexPolygon is convex", g::isConvex(P_convex));
	// batchBooleanTest(
	// 	"generatedConvexPolygon does not produce colinear points",
	// 	N,
	// 	[&P_convex, &N](unsigned int n) {
	// 		return !g::isColinear(
	// 			P_convex[n > 0 ? n - 1 : N - 1], P_convex[n], P_convex[(n + 1) % N]
	// 		);
	// 	}
	// );
	/*
	Test that convexity holds for both clockwise and anticlockwise oriented polygons.
	*/
	// T::Polygon square_clockwise(4);
	// T::Polygon square_anti(4);
	// square_clockwise[0], square_anti[0] = T::Point(0., 0.);
	// square_clockwise[1], square_anti[3] = T::Point(0., 1.);
	// square_clockwise[2], square_anti[2] = T::Point(1., 1.);
	// square_clockwise[3], square_anti[1] = T::Point(1., 0.);
	// booleanTest("isConvex holds for counter-clockwise ordered polygons",
	// g::isConvex(square_anti)); booleanTest("isConvex holds for clockwise ordered polygons",
	// g::isConvex(square_clockwise));
	/*
	Test normaliseConvexPolygon.
	*/
	// booleanTest(
	// 	"normaliseConvexPolygon produces a polygon on the unit interval.",
	// 	g::largestVector(g::normaliseConvexPolygon(P_convex)).first == 1.
	// );
	return 0;
}
