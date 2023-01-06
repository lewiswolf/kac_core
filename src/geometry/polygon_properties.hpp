/*
Utility functions for working with polygons.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <math.h>
#include <utility>

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	T::Point centroid(const T::Polygon& P, const double& area) {
		/*
		This algorithm is used to calculate the geometric centroid of a 2D
		polygon. See http://paulbourke.net/geometry/polygonmesh/ 'Calculating
		the area and centroid of a polygon'.

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

	bool isColinear(const T::Point& a, const T::Point& b, const T::Point& c) {
		/*
		Determines whether or not a given set of three vertices are colinear.
		*/

		return (c.y - b.y) * (b.x - a.x) == (b.y - a.y) * (c.x - b.x);
	}

	bool isConvex(const T::Polygon& P) {
		/*
		Tests whether or not a given array of vertices forms a convex polygon.
		This is achieved using the resultant sign of the cross product for each
		vertex:
			[(x_i - x_i-1), (y_i - y_i-1)] x [(x_i+1 - x_i), (y_i+1 - y_i)]
		See => http://paulbourke.net/geometry/polygonmesh/ 'Determining whether
		or not a polygon (2D) has its vertices ordered clockwise or
		counter-clockwise'.
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

	bool isPointInsidePolygon(const T::Point& p, const T::Polygon& P) {
		/*
		Determines whether or not a cartesian pair is within a polygon.
		Adapted from collidePointPoly() => https://github.com/bmoren/p5.collide2D
		*/

		bool collision = false;
		const unsigned long N = P.size();
		// go through each of the vertices, plus the next vertex in the list
		for (unsigned long n = 0; n < N; n++) {
			T::Point a = P[n];
			T::Point b = P[(n + 1) % N];
			// compare position, flip 'collision' variable back and forth
			if (((a.y >= p.y && b.y < p.y) || (a.y < p.y && b.y >= p.y))
				&& p.x < ((b.x - a.x) * (p.y - a.y)) / (b.y - a.y) + a.x) {
				collision = !collision;
			}
		}
		return collision;
	}

	std::pair<double, std::pair<int, int>> largestVector(const T::Polygon& P) {
		/*
		This function tests each pair of vertices in a given polygon to find the
		largest vector, and returns the length of the vector and its indices.
		*/

		const unsigned long N = P.size();
		double vec_max = 0.;
		long index_i = 0;
		long index_j = 0;
		for (unsigned long i = 0; i < N; i++) {
			for (unsigned long j = i + 1; j < N; j++) {
				double vec = sqrt(pow(P[i].x - P[j].x, 2) + pow(P[i].y - P[j].y, 2));
				if (vec > vec_max) {
					vec_max = vec;
					index_i = i;
					index_j = j;
				}
			}
		}
		return std::make_pair(vec_max, std::make_pair(index_i, index_j));
	}

	double polygonArea(const T::Polygon& P) {
		/*
		An implementation of the shoelace algorithm, first described by Albrecht
		Ludwig Friedrich Meister, which is used to calculate the area of a
		polygon.
		*/

		const unsigned long N = P.size();
		double out = 0.;
		for (unsigned long n = 0; n < N; n++) {
			out += P[n].x * P[(n + 1) % N].y - P[n].y * P[(n + 1) % N].x;
		}
		return abs(out) / 2.;
	}

}