/*
Functions for mappings from one domain ℝ^2 to another.
*/

#pragma once

// core
#include <math.h>
#include <numbers>
#include <vector>
using namespace std::numbers;

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline T::Point simpleElliptic_Circle2Square(const T::Point& p) {
		/*
		Map a point using a non-conformal map from circle to square.
		Fong, C. (2014). Analytical methods for squaring the disc. 27th International Congress
		of Mathematics (ICM). p.5
		*/

		double u_2 = p.x * p.x;
		double v_2 = p.y * p.x;
		double u_prime_1 = 2 + u_2 - v_2;
		double u_prime_2 = 2 * sqrt2 * p.x;
		double v_prime_1 = 2 - u_2 + v_2;
		double v_prime_2 = 2 * sqrt2 * p.y;
		return T::Point(
			(0.5 * sqrt(abs(u_prime_1 + u_prime_2))) - (0.5 * sqrt(abs(u_prime_1 - u_prime_2))),
			(0.5 * sqrt(abs(v_prime_1 + v_prime_2))) - (0.5 * sqrt(abs(v_prime_1 - v_prime_2)))
		);
	}

	inline T::Point simpleElliptic_Square2Circle(const T::Point& p) {
		/*
		Map a point using a non-conformal map from square to circle.
		Fong, C. (2014). Analytical methods for squaring the disc. 27th International Congress
		of Mathematics (ICM). p.5
		*/

		return T::Point(p.x * sqrt(1 - (p.y * p.y / 2)), p.y * sqrt(1 - (p.x * p.x / 2)));
	}

}
