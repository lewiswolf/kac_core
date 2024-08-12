/*
Tests for /geometry.
*/

// core
#include <numbers>
using namespace std::numbers;

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
			return (P_convex[n].x - P_copy[n].x) == 0. && (P_convex[n].y - P_copy[n].y) == 0.;
		}
	);

	/*
	Test the properties of generateConvexPolygon.
	*/
	booleanTest("generatedConvexPolygon produces n vertices", P_convex.size() == N);
	booleanTest("generatedConvexPolygon is convex", g::isConvex(P_convex));
	// batchBooleanTest(
	// 	"generatedConvexPolygon does not produce colinear points",
	// 	N,
	// 	[&P_convex](const unsigned long& n) {
	// 		return !g::isColinear(P_convex[n > 0 ? n - 1 : 9], P_convex[n], P_convex[(n + 1) %
	// 10]);
	// 	}
	// );

	/*
	Test that polygon properties holds for both clockwise and anticlockwise oriented polygons.
	*/
	T::Polygon square_clockwise = {
		T::Point(0., 0.), T::Point(0., 1.), T::Point(1., 1.), T::Point(1., 0.)
	};
	T::Polygon square_anti = {
		T::Point(0., 0.), T::Point(1., 0.), T::Point(1., 1.), T::Point(0., 1.)
	};
	booleanTest("isConvex holds for counter-clockwise ordered polygons", g::isConvex(square_anti));
	booleanTest("isConvex holds for clockwise ordered polygons", g::isConvex(square_clockwise));
	booleanTest("largestVector works anticlockwise", g::largestVector(square_anti).first == sqrt2);
	booleanTest("largestVector works clockwise", g::largestVector(square_clockwise).first == sqrt2);

	/*
	Test that isPointInsideConvexPolygon and isPointInsidePolygon are accurate.
	*/
	booleanTest(
		"isPointInsideConvexPolygon holds.",
		g::isPointInsideConvexPolygon(T::Point(0.5, 0.5), square_clockwise)
	);
	booleanTest(
		"isPointInsidePolygon holds.", g::isPointInsidePolygon(T::Point(0.5, 0.5), square_clockwise)
	);

	/*
	Test that _polygonCentroid works for negative values.
	*/
	T::Polygon square_negative = {
		T::Point(-11., -10.), T::Point(-10., -9.), T::Point(-9., -10.), T::Point(-10., -11.)
	};
	T::Point centroid = g::polygonCentroid(square_negative);
	booleanTest(
		"polygonCentroid holds for negative centroids", centroid.x == -10. && centroid.y == -10.
	);

	/*
	Test normaliseConvexPolygon.
	*/
	// booleanTest(
	// 	"normaliseConvexPolygon produces a polygon on the unit interval.",
	// 	g::largestVector(g::normaliseConvexPolygon(P_convex)).first == 1.
	// );

	/*
	Test Encyclopedia of Triangle Centers.
	*/
	T::Polygon tri = {T::Point(0., 0.), T::Point(1., 0.), T::Point(1., 1.)};
	T::Point incenter = g::ETC::incenter(tri);
	booleanTest(
		"X(1) produces the correct output.",
		(incenter.x - 0.707107) < 0.001 && (incenter.y - 0.292893) < 0.001
	);
	centroid = g::ETC::centroid(tri);
	booleanTest(
		"X(2) produces the correct output.",
		(centroid.x - (2. / 3.)) < 0.001 && (centroid.y - (1. / 3.)) < 0.001
	);
	T::Point circumcenter = g::ETC::circumcenter(tri);
	booleanTest(
		"X(3) produces the correct output.", (circumcenter.x == 0.5) && (circumcenter.y == 0.5)
	);
	T::Point orthocenter = g::ETC::orthocenter(tri);
	booleanTest(
		"X(4) produces the correct output.", (orthocenter.x == 1.) && (orthocenter.y == 0.)
	);

	/*
	Test isPointOnLine is accurate.
	*/
	T::Line test_line = T::Line(T::Point(0., 0.), T::Point(1., 1.));
	booleanTest(
		"isPointOnLine identifies a point on the line.",
		g::isPointOnLine(T::Point(0.5, 0.5), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point on the line.",
		g::isPointOnLine(T::Point(0., 0.), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point on the line.",
		g::isPointOnLine(T::Point(1., 1.), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point not on the line.",
		!g::isPointOnLine(T::Point(0.501, 0.5), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point not on the line.",
		!g::isPointOnLine(T::Point(0.5, 0.501), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point not on the line.",
		!g::isPointOnLine(T::Point(1.001, 1.001), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point not on the line.",
		!g::isPointOnLine(T::Point(-0.001, -0.001), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point not on the line.",
		!g::isPointOnLine(T::Point(-1., -1.), test_line)
	);
	booleanTest(
		"isPointOnLine identifies a point not on the line.",
		!g::isPointOnLine(T::Point(2., 2.), test_line)
	);

	return 0;
}
