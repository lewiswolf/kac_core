/*
Functions for mappings from one domain ℝ^2 to another.
*/

#pragma once

// core
#include <cmath>
#include <numbers>

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	enum class SquareToCircleMethod { Elliptic };

	inline T::Point circleToSquare(
		const T::Point& p, SquareToCircleMethod method = SquareToCircleMethod::Elliptic
	) {
		/*
		Map a point using a non-conformal map from circle to square.
		Fong, C. (2014). Analytical methods for squaring the disc. 27th International Congress
		of Mathematics (ICM). p.5
		*/

		if (method == SquareToCircleMethod::Elliptic) {
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
		return T::Point(0., 0.);
	}

	inline T::Point squareToCircle(
		const T::Point& p, SquareToCircleMethod method = SquareToCircleMethod::Elliptic
	) {
		/*
		Map a point using a non-conformal map from square to circle.
		Fong, C. (2014). Analytical methods for squaring the disc. 27th International Congress
		of Mathematics (ICM). p.5
		*/

		if (method == SquareToCircleMethod::Elliptic) {
			return T::Point(
				p.x * std::sqrt(1 - (p.y * p.y * 0.5)), p.y * std::sqrt(1. - (p.x * p.x * 0.5))
			);
		}
		return T::Point(0., 0.);
	}

	enum class SquareToTriangleMethod { Heitz };

	inline T::Point squareToTriangle(
		const T::Point& p, SquareToTriangleMethod method = SquareToTriangleMethod::Heitz
	) {
		/*
		Map a point on a square to a right angled triangle.
		See:
			Heitz, E. (2019). A low-distortion map between triangle and square
			https://hal.science/hal-02073696v2
		*/

		if (method == SquareToTriangleMethod::Heitz) {
			if (p.y > p.x) {
				return T::Point(p.x * 0.5, p.y - (p.x * 0.5));
			} else {
				return T::Point(p.x - (p.y * 0.5), p.y * 0.5);
			}
		}
		return T::Point(0., 0.);
	}

	inline T::Point triangleToSquare(
		const T::Point& p, SquareToTriangleMethod method = SquareToTriangleMethod::Heitz
	) {
		/*
		Inversion formulas for squareToTriangle.
		*/

		if (method == SquareToTriangleMethod::Heitz) {
			if (p.y > p.x) {
				return T::Point(p.x * 2., p.y + p.x);
			} else {
				return T::Point(p.x + p.y, p.y * 2.);
			}
		}
		return T::Point(0., 0.);
	}

}
