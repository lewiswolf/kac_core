/*
Utility functions for working with polygons.
*/

#pragma once

// core
#include <algorithm>
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
			[(x_i - x_i-1), (y_i - y_i-1)] × [(x_i+1 - x_i), (y_i+1 - y_i)]
		See => http://paulbourke.net/geometry/polygonmesh/ 'Determining whether or not a polygon
		(2D) has its vertices ordered clockwise or counter-clockwise'.
		*/

		// cross product - z component only, see np.cross =>
		// https://numpy.org/doc/stable/reference/generated/numpy.cross.html
		auto crossProductZ = [](T::Point p, T::Point p_plus, T::Point p_minus) {
			return (p.x - p_minus.x) * (p_plus.y - p.y) - (p_plus.x - p.x) * (p.y - p_minus.y);
		};
		// determine the direction of the initial point using the cross product
		const unsigned long N = P.size();
		bool clockwise = crossProductZ(P[0], P[1], P[N - 1]) < 0;
		// loop over remaining points
		for (unsigned long n = 1; n < N; n++) {
			if (crossProductZ(P[n], P[(n + 1) % N], P[n - 1]) < 0 != clockwise) {
				return false;
			}
		}
		return true;
	}

	inline bool isPointInsideConvexPolygon(const T::Point& p, T::Polygon P) {
		/*
		Determines whether or not a cartesian pair is within a polygon, including boundaries.
		Solution 3 => http://paulbourke.net/geometry/polygonmesh/
		*/

		// determine if the polygon is ordered clockwise
		const short clockwise =
			(P[1].x - P[0].x) * (P[2].y - P[1].y) - (P[1].y - P[0].y) * (P[2].x - P[1].x) > 0 ? -1
																							  : 1;
		// go through each of the vertices, plus the next vertex in the list
		const unsigned long N = P.size();
		for (unsigned long n = 0; n < N; n++) {
			T::Point a = P[n];
			T::Point b = P[(n + 1) % N];
			if (((b.x - a.x) * (p.y - b.y) - (b.y - a.y) * (p.x - b.x)) * clockwise > 0.) {
				return false;
			}
		}
		return true;
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
				if (intersection_type == "none" || intersection_type == "vertex") {
					continue;
				} else {
					return false;
				}
			}
		}
		return true;
	}

	inline std::pair<double, std::pair<unsigned long, unsigned long>>
	largestVector(const T::Polygon& P) {
		/*
		This function tests each pair of vertices in a given polygon to find the largest vector, and
		returns the length of the vector and its indices.
		*/

		const unsigned long N = P.size();
		unsigned long index_i = 0;
		unsigned long index_j = 0;
		double vec_max = 0.;
		for (unsigned long i = 0; i < N; i++) {
			for (unsigned long j = i + 1; j < N; j++) {
				double vec = sqrt(pow(P[i].x - P[j].x, 2) + pow(P[i].y - P[j].y, 2));
				if (vec > vec_max) {
					index_i = i;
					index_j = j;
					vec_max = vec;
				}
			}
		}
		return std::make_pair(vec_max, std::make_pair(index_i, index_j));
	}

	inline double polygonArea(const T::Polygon& P) {
		/*
		An implementation of the polygon area algorithm derived using Green's Theorem.
		https://math.blogoverflow.com/2014/06/04/greens-theorem-and-area-of-polygons/
		*/

		const unsigned long N = P.size();
		double out = 0.;
		for (unsigned long n = 0; n < N; n++) {
			out += (P[(n + 1) % N].x + P[n].x) * (P[(n + 1) % N].y - P[n].y);
		}
		return abs(out) / 2.;
	}

	inline T::Point polygonCentroid(const T::Polygon& P, const double& area) {
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
			(x, y) = (
				1/6A * Σ (x_n + x_n+1)(x_n * y_n+1 - x_n+1 * y_n),
				1/6A * Σ (y_n + y_n+1)(x_n * y_n+1 - x_n+1 * y_n),
			)
		 */

		const unsigned long N = P.size();
		double out_x = 0.;
		double out_y = 0.;
		if (N == 3) {
			// Triangles have a much simpler formula, and so these are
			// calculated separately.
			for (unsigned long n = 0; n < 3; n++) {
				out_x += P[n].x;
				out_y += P[n].y;
			}
			return T::Point(out_x / 3., out_y / 3.);
		}
		for (unsigned long n = 0; n < N; n++) {
			double out = (P[n].x * P[(n + 1) % N].y - P[(n + 1) % N].x * P[n].y);
			out_x += (P[n].x + P[(n + 1) % N].x) * out;
			out_y += (P[n].y + P[(n + 1) % N].y) * out;
		}
		return T::Point(abs(out_x) / (6 * area), abs(out_y) / (6 * area));
	}

}
