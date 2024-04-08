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

	// T::Point centroid(const T::Polygon& P) {
	// 	/*
	// 	X(2) Centroid
	// 	*/

	// 	return T::Point();
	// }

}
