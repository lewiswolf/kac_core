/*
Utility functions for working with lines and curves.
*/

#pragma once

// core
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline bool isColinear(const T::Point& a, const T::Point& b, const T::Point& c) {
		/*
		Determines whether or not a given set of three vertices are colinear.
		*/

		return (b.x - a.x) * (c.y - b.y) == (b.y - a.y) * (c.x - b.x);
	}

	inline bool isPointOnLine(const T::Point& p, const T::Line& A) {
		/*
		Determines whether or not a point lies on a line segment.
		*/

		// check if p is within the bounding box of A
		if (p.x < std::min(A.a.x, A.b.x) || p.x > std::max(A.a.x, A.b.x)
			|| p.y < std::min(A.a.y, A.b.y) || p.y > std::max(A.a.y, A.b.y)) {
			return false;
		}
		// check if p is collinear with A.a and A.b
		return isColinear(A.a, A.b, p);
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
		if (isColinear(A.a, A.b, B.a) && isColinear(A.a, A.b, B.b)) {
			if (isPointOnLine(A.a, B) || isPointOnLine(A.b, B) || isPointOnLine(B.a, A)
				|| isPointOnLine(B.b, A)) {
				return std::make_pair(
					"colinear",
					T::Point(
						(A.a.x + A.b.x + B.a.x + B.b.x) / 4, (A.a.y + A.b.y + B.a.y + B.b.y) / 4
					)
				);
			}
		} else {
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
		}
		// return the null case
		return std::make_pair("none", T::Point());
	}

	inline T::Point lineMidpoint(const T::Line& L) {
		/*
		Find the midpoint of a line.
		*/

		return T::Point((L.a.x + L.b.x) * 0.5, (L.a.y + L.b.y) * 0.5);
	}

}
