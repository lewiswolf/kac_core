/*
Utility functions for working with lines and curves.
*/

#pragma once

// core
#include <string>
#include <utility>

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline T::BooleanImage bresenham(T::BooleanImage& M, const T::Line& L) {
		/*
		Apply the Bresenham line drawing algorithm to an input matrix.
		input:
			M = input matrix.
			L = line to draw, such that x ∈ [0, 1] && y ∈ [0, 1].
		*/

		// assert line is within the unit interval
		if (L.a.x > 1. || L.a.x < 0. || L.a.y > 1. || L.a.y < 0. || L.b.x > 1. || L.b.x < 0.
			|| L.b.y > 1. || L.b.y < 0.) {
			throw std::invalid_argument(
				"The line L must be within the unit interval, such that x ∈ [0, 1] && y ∈ [0, 1]."
			);
		}
		// handle discretisation
		const unsigned long x_0 = lround(L.a.x * (M.size() - 1));
		const unsigned long y_0 = lround(L.a.y * (M[0].size() - 1));
		long dx = lround(L.b.x * (M.size() - 1)) - x_0;
		long dy = lround(L.b.y * (M[0].size() - 1)) - y_0;
		// configure directions
		short xx, xy, yx, yy;
		if (std::abs(dx) > std::abs(dy)) {
			xx = dx > 0 ? 1 : -1;
			xy = 0;
			yx = 0;
			yy = dy > 0 ? 1 : -1;
		} else {
			std::swap(dx, dy);
			xx = 0;
			xy = dx > 0 ? 1 : -1;
			yx = dy > 0 ? 1 : -1;
			yy = 0;
		}
		dx = std::abs(dx);
		dy = std::abs(dy);
		// paint line
		unsigned long y = 0;
		long D = 2 * dy - dx;
		for (unsigned long x = 0; x < (dx + 1); x++) {
			M[x_0 + x * xx + y * yx][y_0 + x * xy + y * yy] = 1;
			// reposition y
			if (D >= 0) {
				y += 1;
				D -= 2 * dx;
			}
			D += 2 * dy;
		}
		return M;
	}

	inline bool isColinear(const T::Point& a, const T::Point& b, const T::Point& c) {
		/*
		Determines whether or not a given set of three vertices are colinear.
		*/

		return (c.y - b.y) * (b.x - a.x) == (b.y - a.y) * (c.x - b.x);
	}

	inline std::pair<std::string, T::Point> lineIntersection(const T::Line& A, const T::Line& B) {
		/*
		This function determines whether a line has an intersection, and returns it's type as well
		as the point of intersection (if one exists).
		input:
			A, B - Line segments to compare.
		output:
			type -
				'none'		No intersection.
				'intersect' The general case where lines intersect one another.
				'vertex'	This is the special case when two lines share a vertex.
				'branch'	This is the special case when a vertex lies within another line. For
							example, B creates an intersection at point B.a when B.a lies on the
							open interval (A.a, A.b).
				'colinear'	This is the special case when the two lines overlap.
			point -
				'none'		Empty point.
				'intersect' The point of intersection ∈ (A.a, A.b) & (B.a, B.b).
				'vertex'	The shared vertex.
				'branch'	The branching vertex.
				'colinear'	The midpoint between all 4 vertices.
		*/

		// search for shared vertices
		if (A.a.x == B.a.x && A.a.y == B.a.y || A.a.x == B.b.x && A.a.y == B.b.y) {
			return std::make_pair("vertex", A.a);
		} else if (A.b.x == B.a.x && A.b.y == B.a.y || A.b.x == B.b.x && A.b.y == B.b.y) {
			return std::make_pair("vertex", A.b);
		}
		// test for colinear cases.
		unsigned short colinearities = 0;
		colinearities += isColinear(A.a, A.b, B.a);
		colinearities += isColinear(A.b, B.a, B.b);
		colinearities += isColinear(B.a, B.b, A.a);
		colinearities += isColinear(B.b, A.a, A.b);
		if (colinearities == 4) {
			return std::make_pair(
				"colinear",
				T::Point((A.a.x + A.b.x + B.a.x + B.b.x) / 4, (A.a.y + A.b.y + B.a.y + B.b.y) / 4)
			);
		}
		// calculate the general case using distance to intersection point.
		double u_A = ((B.b.x - B.a.x) * (A.a.y - B.a.y) - (B.b.y - B.a.y) * (A.a.x - B.a.x))
				   / ((B.b.y - B.a.y) * (A.b.x - A.a.x) - (B.b.x - B.a.x) * (A.b.y - A.a.y));
		double u_B = ((A.b.x - A.a.x) * (A.a.y - B.a.y) - (A.b.y - A.a.y) * (A.a.x - B.a.x))
				   / ((B.b.y - B.a.y) * (A.b.x - A.a.x) - (B.b.x - B.a.x) * (A.b.y - A.a.y));
		if (u_A >= 0 && u_A <= 1 && u_B >= 0 && u_B <= 1) {
			T::Point p = T::Point(A.a.x + u_A * (A.b.x - A.a.x), A.a.y + u_A * (A.b.y - A.a.y));
			// test for adjacent case
			if (A.a.x == p.x && A.a.y == p.y) {
				return std::make_pair("adjacent", A.a);
			} else if (A.b.x == p.x && A.b.y == p.y) {
				return std::make_pair("adjacent", A.b);
			} else if (B.a.x == p.x && B.a.y == p.y) {
				return std::make_pair("adjacent", B.a);
			} else if (B.b.x == p.x && B.b.y == p.y) {
				return std::make_pair("adjacent", B.b);
			}
			// return general case
			return std::make_pair("intersect", p);
		}
		// return the null case
		return std::make_pair("none", T::Point());
	}

	inline T::Point lineMidpoint(const T::Line& L) {
		/*
		Find the midpoint of a line.
		*/

		return T::Point((L.a.x + L.b.x) / 2., (L.a.y + L.b.y) / 2.);
	}

}
