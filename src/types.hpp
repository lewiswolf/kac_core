/*
🙇‍♀️🙇🙇‍♂️
*/

#pragma once

// core
#include <array>
#include <math.h>
#include <utility>
#include <vector>

namespace kac_core::types {

	// Matrices
	typedef std::vector<double> Matrix_1D;
	typedef std::vector<std::vector<double>> Matrix_2D;
	typedef std::vector<std::vector<short>> BooleanImage;

	typedef struct Point {
		/*
		A point on the Euclidean plane.
		*/

		// vars
		double x = 0.;
		double y = 0.;
		double r() { return hypot(x, y); }
		double theta() { return atan2(y, x); }

		// constructors
		Point() {}
		Point(double x, double y): x(x), y(y) {}
		Point(std::array<double, 2> _array): x(_array[0]), y(_array[1]) {}
		Point(std::pair<double, double> _pair): x(_pair.first), y(_pair.second) {}

		// methods
		void updatePol(double r, double theta) {
			/*
			Update the point using polar coordinates.
			*/

			x = r * cos(theta);
			y = r * sin(theta);
		}
	} Point;

	typedef struct Line {
		/*
		A straight line from point a to point b.
		*/

		// vars
		Point a;
		Point b;

		// constructors
		Line() {}
		Line(Point a, Point b): a(a), b(b) {}
	} Line;

	// A polygon defined on the Euclidean plane.
	typedef std::vector<Point> Polygon;

}
