/*
Various formulas according to the Encyclopedia of Triangle Centers.
https://faculty.evansville.edu/ck6/encyclopedia/ETC.html
*/

#pragma once

// core
#include <stdexcept>

// src
#include "../types.hpp"
namespace T = kac_core::types;

void typeGuard(const T::Polygon& P) {
	/*
	Enforce that each function is performed on a triangle.
	*/
	if (P.size() != 3) {
		throw std::invalid_argument(
			"Triangle centers can only be calculated for three sided polygons."
		);
	}
}

namespace kac_core::geometry::ETC {

	T::Point incenter(const T::Polygon& P) {
		/*
		X(1) Incenter
		output:
			( x, y ) = coordinates of the incenter.
		*/

		typeGuard(P);
		double a = sqrt(pow(P[1].x - P[2].x, 2) + pow(P[1].y - P[2].y, 2));
		double b = sqrt(pow(P[0].x - P[2].x, 2) + pow(P[0].y - P[2].y, 2));
		double c = sqrt(pow(P[0].x - P[1].x, 2) + pow(P[0].y - P[1].y, 2));
		return T::Point(
			(a * P[0].x + b * P[1].x + c * P[2].x) / (a + b + c),
			(a * P[0].y + b * P[1].y + c * P[2].y) / (a + b + c)
		);
	}

	T::Point centroid(const T::Polygon& P) {
		/*
		X(2) Centroid
		output:
			( x, y ) = coordinates of the centroid.
		*/

		typeGuard(P);
		return T::Point((P[0].x + P[1].x + P[2].x) / 3., (P[0].y + P[1].y + P[2].y) / 3.);
	}

	T::Point circumcenter(const T::Polygon& P) {
		/*
		X(2) Circumcenter
		output:
			( x, y ) = coordinates of the circumcenter.
		*/

		typeGuard(P);
		const double a_2 = P[0].x * P[0].x + P[0].y * P[0].y;
		const double b_2 = P[1].x * P[1].x + P[1].y * P[1].y;
		const double c_2 = P[2].x * P[2].x + P[2].y * P[2].y;
		const double d = 2
					   * (P[0].x * (P[1].y - P[2].y) + P[1].x * (P[2].y - P[0].y)
						  + P[2].x * (P[0].y + P[1].y));
		return T::Point(
			((a_2 * (P[1].y - P[2].y)) + (b_2 * (P[2].y - P[0].y)) + (c_2 * (P[0].y - P[1].y))) / d,
			((a_2 * (P[2].x - P[1].x)) + (b_2 * (P[0].x - P[2].x)) + (c_2 * (P[1].x - P[0].x))) / d
		);
	}

	T::Point orthocenter(const T::Polygon& P) {
		/*
		X(4) Orthocenter
		output:
			( x, y ) = coordinates of the orthocenter.
		*/

		typeGuard(P);
		const double a = P[1].x * (P[0].x - P[2].x) + P[1].y * (P[0].y - P[2].y);
		const double b = P[0].x * (P[1].x - P[2].x) + P[0].y * (P[1].y - P[2].y);
		const double c = (P[2].x - P[1].x) * (P[2].y - P[0].y);
		const double d = (P[2].y - P[1].y) * (P[2].x - P[0].x);
		return T::Point(
			(a * (P[2].y - P[1].y) - b * (P[2].y - P[0].y)) / (c - d),
			(a * (P[2].x - P[1].x) - b * (P[2].x - P[0].x)) / (d - c)
		);
	}

}
