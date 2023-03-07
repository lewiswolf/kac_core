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
	int N = 10;
	T::Polygon p = g::generateConvexPolygon(N);
	booleanTest("generatedConvexPolygon produces n vertices", p.size() == N);
	booleanTest(
		"convexNormalisation produces a polygon on the unit interval.",
		g::largestVector(g::convexNormalisation(p)).first == 1.
	);
	booleanTest("generatedConvexPolygon is convex", g::isConvex(p));
	batchBooleanTest(
		"generatedConvexPolygon does not produce colinear points",
		N,
		[&N, &p](unsigned int i) {
			return !g::isColinear(p[i > 0 ? i - 1 : N - 1], p[i], p[(i + 1) % N]);
		}
	);

	T::Polygon seed_test = g::generateConvexPolygon(10, 1);
	for (unsigned int i = 0; i < 10; i++) {
		std::cout << " " << seed_test[i].x << " " << seed_test[i].y << std::endl;
	}
	T::Polygon seed_expected(10);
	T::Matrix_2D seed_m = {
		{0.470347, 0.459701},
		{-0.166557, 0.427046},
		{-0.359816, 0.268253},
		{-0.491853, -0.357375},
		{-0.430073, -0.459701},
		{-0.0772344, -0.423636},
		{0.0708572, -0.362168},
		{0.373845, -0.166994},
		{0.478783, -0.0225504},
		{0.491853, 0.159062}};
	for (unsigned int i = 0; i < 10; i++) {
		seed_expected[i] = T::Point(seed_m[i][0], seed_m[i][1]);
	}
	batchBooleanTest(
		"generatedConvexPolygon can be correctly controlled using the seed "
		"parameter",
		10,
		[&seed_test, &seed_expected](unsigned int i) {
			return abs(seed_test[i].x - seed_expected[i].x) < 0.01
				&& abs(seed_test[i].y - seed_expected[i].y) < 0.01;
		}
	);

	T::Polygon square_clockwise(4);
	T::Polygon square_anti(4);
	square_clockwise[0], square_anti[0] = T::Point(0.0, 0.0);
	square_clockwise[1], square_anti[3] = T::Point(0.0, 1.0);
	square_clockwise[2], square_anti[2] = T::Point(1.0, 1.0);
	square_clockwise[3], square_anti[1] = T::Point(1.0, 0.0);
	booleanTest("isConvex works for counter-clockwise ordered polygons", g::isConvex(square_anti));
	booleanTest("isConvex works for clockwise ordered polygons", g::isConvex(square_clockwise));

	return 0;
}