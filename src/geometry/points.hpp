/*
Functions for manipulating or mapping points on the Euclidean plane.
*/

#pragma once

// core
#include <cmath>
#include <stdexcept>

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline T::Point rotatePoint(const T::Point& p, const double& theta) {
		const double cos_theta = std::cos(theta);
		const double sin_theta = std::sin(theta);
		return T::Point(
			(p.x * cos_theta) - (p.y * sin_theta), (p.x * sin_theta) + (p.y * cos_theta)
		);
	}

	inline std::array<double, 2> cartesianToPolar(const T::Point& p) {
		return {std::hypot(p.x, p.y), std::atan2(p.y, p.x)};
	}

	inline T::Point polarToCartesian(const double& r, const double& theta) {
		return T::Point(r * std::cos(theta), r * std::sin(theta));
	}

	inline std::array<double, 3> cartesianToTrilinear(const T::Point& p, const T::Polygon& P) {
		/*
		Convert a cartesian coordinate to a trilinear coordinate relative to a given triangle.
		*/

		if (P.size() != 3) {
			throw std::invalid_argument(
				"Trilinear coordinates can only be calculated for three sided polygons."
			);
		}
		const auto pointToLineDistance = [p](T::Point a, T::Point b) -> double {
			const double A = a.y - b.y;
			const double B = b.x - a.x;
			const double C = a.x * b.y - b.x * a.y;
			return std::abs(A * p.x + B * p.y + C) / std::hypot(A, B);
		};
		const double alpha = pointToLineDistance(P[1], P[2]);
		const double beta = pointToLineDistance(P[2], P[0]);
		const double gamma = pointToLineDistance(P[0], P[1]);
		return {alpha, beta, gamma};
	}

	T::Point
	trilinearToCartesian(const double& u, const double& v, const double& w, const T::Polygon& P) {
		/*
		Convert a trilinear coordinate to a cartesian coordinate relative to a given triangle.
		*/

		if (P.size() != 3) {
			throw std::invalid_argument(
				"Trilinear coordinates can only be calculated for three sided polygons."
			);
		}
		// compute side lengths opposite vertices
		const double a = std::hypot(P[2].x - P[1].x, P[2].y - P[1].y);
		const double b = std::hypot(P[0].x - P[2].x, P[0].y - P[2].y);
		const double c = std::hypot(P[1].x - P[0].x, P[1].y - P[0].y);
		// convert to barycentric coordinates
		const double sum = 1. / (a * u + b * v + c * w);
		const double lambda_1 = a * u * sum;
		const double lambda_2 = b * v * sum;
		const double lambda_3 = c * w * sum;
		// convert to cartesian coordinates
		return T::Point(
			lambda_1 * P[0].x + lambda_2 * P[1].x + lambda_3 * P[2].x,
			lambda_1 * P[0].y + lambda_2 * P[1].y + lambda_3 * P[2].y
		);
	}

}
