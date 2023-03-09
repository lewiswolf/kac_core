/*
Utility functions for working with lines and curves.
*/

#pragma once

// core
#include <utility>

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline std::pair<bool, T::Point> lineIntersection(const T::Line& A, const T::Line& B) {
		/*
		Finds the point at which two lines intersect.
		collisionLineLine() => https://github.com/bmoren/p5.collide2D
		*/
		// calculate the distance to intersection point
		double u_A = ((B.b.x - B.a.x) * (A.a.y - B.a.y) - (B.b.y - B.a.y) * (A.a.x - B.a.x))
				   / ((B.b.y - B.a.y) * (A.b.x - A.a.x) - (B.b.x - B.a.x) * (A.b.y - A.a.y));
		double u_B = ((A.b.x - A.a.x) * (A.a.y - B.a.y) - (A.b.y - A.a.y) * (A.a.x - B.a.x))
				   / ((B.b.y - B.a.y) * (A.b.x - A.a.x) - (B.b.x - B.a.x) * (A.b.y - A.a.y));
		// if u_A and u_B are between 0-1, lines are colliding.
		if (u_A >= 0 && u_A <= 1 && u_B >= 0 && u_B <= 1) {
			return std::make_pair(
				true, T::Point(A.a.x + u_A * (A.b.x - A.a.x), A.a.y + u_A * (A.b.y - A.a.y))
			);
		} else {
			return std::make_pair(false, T::Point());
		}
	}

}