/*
Functions for mappings from one domain ℝ^2 to another.
*/

#pragma once

// core
#include <cmath>
#include <numbers>
#include <vector>

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

		const double u_2 = p.x * p.x;
		const double v_2 = p.y * p.y;
		const double u_prime_1 = 2. + u_2 - v_2;
		const double u_prime_2 = 2. * std::numbers::sqrt2 * p.x;
		const double v_prime_1 = 2. - u_2 + v_2;
		const double v_prime_2 = 2. * std::numbers::sqrt2 * p.y;
		return T::Point(
			(std::sqrt(std::abs(u_prime_1 + u_prime_2))
			 - std::sqrt(std::abs(u_prime_1 - u_prime_2)))
				* 0.5,
			(std::sqrt(std::abs(v_prime_1 + v_prime_2))
			 - std::sqrt(std::abs(v_prime_1 - v_prime_2)))
				* 0.5
		);
	}

	inline T::Point simpleElliptic_Square2Circle(const T::Point& p) {
		/*
		Map a point using a non-conformal map from square to circle.
		Fong, C. (2014). Analytical methods for squaring the disc. 27th International Congress
		of Mathematics (ICM). p.5
		*/

		return T::Point(
			p.x * std::sqrt(1 - (p.y * p.y * 0.5)), p.y * std::sqrt(1. - (p.x * p.x * 0.5))
		);
	}

}
