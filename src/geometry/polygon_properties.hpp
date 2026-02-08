/*
Utility functions for working with polygons.
*/

#pragma once

// core
#include <algorithm>
#include <cmath>
#include <math.h>
#include <string>
#include <utility>

// src
#include "../types.hpp"
#include "./lines.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline bool isConvex(const T::Polygon& P) {
		/*
		Tests whether or not a given array of vertices forms a convex polygon. This is achieved
		using the resultant sign of the cross product for each vertex:
			(x_n - x_n-1)(y_n+1 - y_n) - (x_n+1 - x_n)(y_n - y_n-1)
		See => http://paulbourke.net/geometry/polygonmesh/ 'Determining whether or not a polygon
		(2D) has its vertices ordered clockwise or counter-clockwise'.
		*/

		// cross product - z component only, see np.cross =>
		// https://numpy.org/doc/stable/reference/generated/numpy.cross.html
		const auto crossProductZ = [](T::Point a, T::Point b, T::Point c) -> double {
			return (b.x - a.x) * (c.y - b.y) - (c.x - b.x) * (b.y - a.y);
		};
		// determine the direction of the initial point using the cross product
		const std::size_t N = P.size();
		const bool clockwise = crossProductZ(P[N - 1], P[0], P[1]) < 0.;
		// loop over remaining points
		for (std::size_t n = 1; n < N; n++) {
			if (crossProductZ(P[n - 1], P[n], P[(n + 1) % N]) < 0. != clockwise) {
				return false;
			}
		}
		return true;
	}

	inline bool isPointInsideConvexPolygon(const T::Point& p, const T::Polygon& P) {
		/*
		Determines whether or not a cartesian pair is within a polygon, including boundaries.
		Solution 3 => http://paulbourke.net/geometry/polygonmesh/
		*/

		auto crossProductZ = [](T::Point a, T::Point b, T::Point p) {
			return (b.x - a.x) * (p.y - a.y) - (p.x - a.x) * (b.y - a.y);
		};
		// determine if the polygon is ordered clockwise
		const short clockwise = crossProductZ(P[0], P[1], P[2]) > 0 ? -1 : 1;
		// go through each of the vertices, and test with p
		const unsigned long N = P.size();
		for (unsigned long n = 0; n < N; n++) {
			if (p.x == P[n].x && p.y == P[n].y) {
				return true;
			}
		}
		// determine if the point is always on the right side of the line
		for (unsigned long n = 0; n < N; n++) {
			if (crossProductZ(P[n], P[(n + 1) % N], p) * clockwise > 0.) {
				return false;
			}
		}
		return true;
	}

	inline bool isPointInsidePolygon(const T::Point& p, const T::Polygon& P) {
		/*
		Determines whether or not a cartesian pair is within a polygon, including boundaries.
		This algorithm builds upon the ray tracing ideas shown in solution 1
			=> https://paulbourke.net/geometry/polygonmesh/
		*/

		const unsigned long N = P.size();
		// create a ray that extends to the right of the polygon
		double max_x = 0.;
		for (unsigned long n = 0; n < N; n++) { max_x = std::max(P[n].x, max_x); }
		T::Line ray = T::Line(p, T::Point(max_x + 1., p.y));
		// count the number of times the ray is intersected
		unsigned long count = 0;
		for (unsigned long n = 0; n < N; n++) {
			// return true if point is a vertex
			if (P[n].x == p.x && P[n].y == p.y) {
				return true;
			}
			// return true if point is on the line
			T::Line A = T::Line(P[n], P[(n + 1) % N]);
			if (isPointOnLine(p, A)) {
				return true;
			}
			// general case
			if (lineIntersection(ray, A).first == "intersect") {
				count++;
			}
		}
		return count % 2 == 1;
	}

	inline bool isSimple(const T::Polygon& P) {
		/*
		Determine if a polygon is simple by checking for intersections.
		*/

		const unsigned long N = P.size();
		for (unsigned long i = 0; i < N - 2; i++) {
			for (unsigned long j = i + 1; j < N; j++) {
				std::string intersection_type =
					lineIntersection(T::Line(P[i], P[i + 1]), T::Line(P[j], P[(j + 1) % N])).first;
				if (intersection_type != "none" && intersection_type != "vertex") {
					return false;
				}
			}
		}
		return true;
	}

	inline std::pair<double, std::pair<std::size_t, std::size_t>>
	largestVector(const T::Polygon& P) {
		/*
		This function tests each pair of vertices in a given polygon to find the largest vector, and
		returns the length of the vector and its indices.
		*/

		const std::size_t N = P.size();
		std::pair<std::size_t, std::size_t> index = {0, 0};
		double vec_max2 = 0.;
		for (std::size_t i = 0; i < N; i++) {
			for (std::size_t j = i + 1; j < N; j++) {
				double dx = P[i].x - P[j].x;
				double dy = P[i].y - P[j].y;
				double vec2 = dx * dx + dy * dy;
				if (vec2 > vec_max2) {
					index.first = i;
					index.second = j;
					vec_max2 = vec2;
				}
			}
		}
		return {std::sqrt(vec_max2), index};
	}

	inline double polygonArea(const T::Polygon& P) {
		/*
		An implementation of the polygon area algorithm derived using Green's Theorem.
		Green's theorem enforces that polygon's with anti-clockwise oriented vertices have positive
		area, and those with clockwise oriented vertices have negative area.
		See: https://math.blogoverflow.com/2014/06/04/greens-theorem-and-area-of-polygons/
		*/

		const std::size_t N = P.size();
		double out = 0.;
		for (std::size_t n = 0; n < N; n++) {
			T::Point P_1 = P[(n + 1) % N];
			out += (P_1.x + P[n].x) * (P_1.y - P[n].y);
		}
		return out * 0.5;
	}

	inline T::Point polygonCentroid(const T::Polygon& P) {
		/*
		This algorithm is used to calculate the geometric centroid of a 2D polygon.
		See http://paulbourke.net/geometry/polygonmesh/ 'Calculating the area and centroid of a
		polygon'.
		output:
			for N == 3 ->
			(x, y) = (
				(x_0 + x_1 + x_2) / 3,
				(y_0 + y_1 + y_2) / 3,
			)
			for N > 3 ->
			A = Σ (x_n * y_n+1) - (x_n+1 * y_n)
			(x, y) = (
				1/3A * Σ (x_n + x_n+1)(x_n * y_n+1 - x_n+1 * y_n),
				1/3A * Σ (y_n + y_n+1)(x_n * y_n+1 - x_n+1 * y_n),
			)
		 */

		const std::size_t N = P.size();
		if (N == 3) {
			// Triangles have a much simpler formula, and so these are
			// calculated separately.
			return T::Point((P[0].x + P[1].x + P[2].x) / 3., (P[0].y + P[1].y + P[2].y) / 3.);
		}
		double area = 0.;
		double out_x = 0.;
		double out_y = 0.;
		for (std::size_t n = 0; n < N; n++) {
			T::Point P_1 = P[(n + 1) % N];
			area += (P_1.x + P[n].x) * (P_1.y - P[n].y);
			double scalar = (P[n].x * P_1.y - P_1.x * P[n].y);
			out_x += (P[n].x + P_1.x) * scalar;
			out_y += (P[n].y + P_1.y) * scalar;
		}
		return T::Point(out_x / (3. * area), out_y / (3. * area));
	}

}
