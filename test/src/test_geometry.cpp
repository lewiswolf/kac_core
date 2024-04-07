/*
Tests for /geometry.
*/

// src
#include <kac_core.hpp>
namespace T = kac_core::types;		 // types
namespace g = kac_core::geometry;	 // geometry

// test
#include "./utils.hpp"

int main() {
	/*
	Initialise polygon.
	*/
	unsigned long N = 10;
	T::Polygon P_convex = g::generateConvexPolygon(N, 1);

	/*
	Test the seed argument of generateConvexPolygon.
	*/
	T::Polygon P_copy = g::generateConvexPolygon(N, 1);
	batchBooleanTest(
		"generatedConvexPolygon produces the expected output for a given seed",
		N,
		[&P_convex, &P_copy](const unsigned long& n) {
			return abs(P_convex[n].x - P_copy[n].x) == 0. && abs(P_convex[n].y - P_copy[n].y) == 0.;
		}
	);

	/*
	Test the properties of generateConvexPolygon.
	*/
	booleanTest("generatedConvexPolygon produces n vertices", P_convex.size() == N);
	booleanTest("generatedConvexPolygon is convex", g::isConvex(P_convex));
	// booleanTest(
	// 	"generatedConvexPolygon does not produce colinear points 9 0 1",
	// 	!g::isColinear(P_convex[9], P_convex[0], P_convex[1])
	// );
	// batchBooleanTest(
	// 	"generatedConvexPolygon does not produce colinear points",
	// 	N,
	// 	[&P_convex](const unsigned long& n) {
	// 		return !g::isColinear(P_convex[n > 0 ? n - 1 : 9], P_convex[n], P_convex[(n + 1) %
	// 10]);
	// 	}
	// );

	/*
	Test that convexity holds for both clockwise and anticlockwise oriented polygons.
	*/
	T::Polygon square_clockwise(4);
	T::Polygon square_anti(4);
	square_clockwise[0], square_anti[0] = T::Point(0., 0.);
	square_clockwise[1], square_anti[3] = T::Point(0., 1.);
	square_clockwise[2], square_anti[2] = T::Point(1., 1.);
	square_clockwise[3], square_anti[1] = T::Point(1., 0.);
	booleanTest("isConvex holds for counter-clockwise ordered polygons", g::isConvex(square_anti));
	booleanTest("isConvex holds for clockwise ordered polygons", g::isConvex(square_clockwise));

	// /*
	// Test normaliseConvexPolygon.
	// */
	// booleanTest(
	// 	"normaliseConvexPolygon produces a polygon on the unit interval.",
	// 	g::largestVector(g::normaliseConvexPolygon(P_convex)).first == 1.
	// );

	return 0;
}
