/*
🙇‍♀️🙇🙇‍♂️
*/

#pragma once

// core
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace kac_core::types {

	// Matrices
	typedef std::vector<double> Matrix_1D;
	typedef std::vector<std::vector<double>> Matrix_2D;
	typedef std::vector<short> BooleanImage_1D;
	typedef std::vector<std::vector<short>> BooleanImage_2D;

	typedef struct Point {
		/*
		A point on the Euclidean plane.
		*/

		// vars
		double x {0.};
		double y {0.};

		// constructors
		constexpr Point() noexcept = default;
		constexpr Point(double _x, double _y) noexcept: x(_x), y(_y) {}
		explicit constexpr Point(const std::array<double, 2>& arr) noexcept: x(arr[0]), y(arr[1]) {}
		explicit constexpr Point(const std::pair<double, double>& p) noexcept
			: x(p.first), y(p.second) {}

		// accessors
		[[nodiscard]] double r() const noexcept { return std::hypot(x, y); }
		[[nodiscard]] double theta() const noexcept { return std::atan2(y, x); }

		// mutators
		void updatePolar(double r, double theta) noexcept {
			/*
			Update the point using polar coordinates.
			*/
			x = r * std::cos(theta);
			y = r * std::sin(theta);
		}
	} Point;

	typedef struct Line {
		/*
		A straight line segment from point a to point b.
		*/

		// vars
		Point a {};
		Point b {};

		// constructors
		constexpr Line() noexcept = default;
		constexpr Line(const Point& _a, const Point& _b) noexcept: a(_a), b(_b) {}
	} Line;

	// A polygon defined on the Euclidean plane.
	typedef std::vector<Point> Polygon;

}
