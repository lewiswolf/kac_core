/*
Functions for mappings from one domain ℝ^2 to another.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <math.h>

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	T::Point simpleElliptic_Circle2Square(const T::Point& p) {
		/*
			Map a point using a non-conformal map from circle to square.
			Fong, C. (2014). Analytical methods for squaring the disc. 27th
			International Congress of Mathematics (ICM). p.5
		*/

		double u_2 = pow(p.x, 2);
		double v_2 = pow(p.y, 2);
		double u_prime = 2 * M_SQRT2 * p.x;
		double v_prime = 2 * M_SQRT2 * p.y;
		return T::Point(
			(0.5 * sqrt(2 + u_2 - v_2 + u_prime)) - (0.5 * sqrt(2 + u_2 - v_2 - u_prime)),
			(0.5 * sqrt(2 - u_2 + v_2 + v_prime)) - (0.5 * sqrt(2 - u_2 + v_2 - v_prime))
		);
	}

	T::Point simpleElliptic_Square2Circle(const T::Point& p) {
		/*
			Map a point using a non-conformal map from square to circle.
			Fong, C. (2014). Analytical methods for squaring the disc. 27th
			International Congress of Mathematics (ICM). p.5
		*/

		return T::Point(p.x * sqrt(1 - (pow(p.y, 2) / 2)), p.y * sqrt(1 - (pow(p.x, 2) / 2)));
	}

}