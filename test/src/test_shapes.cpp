/*
Tests for /shapes.
*/

// src
#include <kac_core.hpp>
using namespace geometry;

// test
#include "./utils.hpp"

int main() {
	int N = 10;
	Vertices v = generateConvexPolygon(N);
	booleanTest("generatedConvexPolygon produces n vertices", v.size() == N);
	booleanTest("generatedConvexPolygon is convex", isConvex(v));
	batchBooleanTest(
		"generatedConvexPolygon does not produce colinear points",
		N,
		[&N, &v](unsigned int i) {
			return !isColinear(v[i > 0 ? i - 1 : N - 1], v[i], v[(i + 1) % N]);
		}
	);

	Vertices seed_test = generateConvexPolygon(10, 1);
	Vertices seed_expected(10);
	Matrix_2D seed_m = {
		{0.288255, 0.411634},
		{-0.244504, 0.0358297},
		{-0.467343, -0.251818},
		{-0.466911, -0.411634},
		{-0.22722, -0.392744},
		{0.0281766, -0.36587},
		{0.115598, -0.35249},
		{0.335813, -0.0359161},
		{0.382849, 0.100085},
		{0.467343, 0.40135}};
	for (unsigned int i = 0; i < 10; i++) {
		seed_expected[i] = Point(seed_m[i][0], seed_m[i][1]);
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

	Vertices p1(4);
	Vertices p2(4);
	p1[0], p2[0] = Point(0.0, 0.0);
	p1[1], p2[3] = Point(0.0, 1.0);
	p1[2], p2[2] = Point(1.0, 1.0);
	p1[3], p2[1] = Point(1.0, 0.0);
	booleanTest(
		"isConvex works for counter-clockwise ordered polygons", isConvex(p1)
	);
	booleanTest("isConvex works for clockwise ordered polygons", isConvex(p2));
	return 0;
}