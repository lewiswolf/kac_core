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
using namespace kac_core::types;

namespace kac_core { namespace geometry {

	Point centroid(const Vertices& V, const double& area) {
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

		const unsigned long N = V.size();
		double out_x = 0.;
		double out_y = 0.;
		if (N == 3) {
			// Triangles have a much simpler formula, and so these are
			// calculated separately.
			for (unsigned long n = 0; n < 3; n++) {
				out_x += V[n].x;
				out_y += V[n].y;
			}
			return Point(out_x / 3., out_y / 3.);
		}
		for (unsigned long n = 0; n < N; n++) {
			double out =
				(V[n].x * V[(n + 1) % N].y - V[(n + 1) % N].x * V[n].y);
			out_x += (V[n].x + V[(n + 1) % N].x) * out;
			out_y += (V[n].y + V[(n + 1) % N].y) * out;
		}
		return Point(abs(out_x) / (6 * area), abs(out_y) / (6 * area));
	}

	bool isColinear(const Point& a, const Point& b, const Point& c) {
		/*
		Determines whether or not a given set of three vertices are
		colinear.
		*/

		return (c.y - b.y) * (b.x - a.x) == (b.y - a.y) * (c.x - b.x);
	}

	bool isConvex(const Vertices& v) {
		/*
		Tests whether or not a given array of vertices forms a convex
		polygon. This is achieved using the resultant sign of the cross
		product for each vertex:
			[(x_i - x_i-1), (y_i - y_i-1)] x [(x_i+1 - x_i), (y_i+1 - y_i)]
		See => http://paulbourke.net/geometry/polygonmesh/ 'Determining
		whether or not a polygon (2D) has its vertices ordered clockwise or
		counter-clockwise'.
		*/

		// cross product - z component only, see np.cross =>
		// https://numpy.org/doc/stable/reference/generated/numpy.cross.html
		auto crossProductZ = [](Point p, Point p_plus, Point p_minus) {
			return (p.x - p_minus.x) * (p_plus.y - p.y)
				- (p_plus.x - p.x) * (p.y - p_minus.y);
		};
		// determine the direction of the initial point using the cross
		// product
		bool clockwise = crossProductZ(v[0], v[1], v[v.size() - 1]) < 0;
		// loop over remaining points
		for (unsigned int i = 1; i < v.size(); i++) {
			if (crossProductZ(v[i], v[(i + 1) % v.size()], v[i - 1]) < 0
				!= clockwise) {
				return false;
			}
		}
		return true;
	}

	std::pair<double, std::pair<int, int>> largestVector(const Vertices& V) {
		/*
		This function tests each pair of vertices in a given polygon to find
		the largest vector, and returns the length of the vector and its
		indices.
		*/

		const unsigned long N = V.size();
		double vec_max = 0.;
		long index_i = 0;
		long index_j = 0;
		for (unsigned long i = 0; i < N; i++) {
			for (unsigned long j = i + 1; j < N; j++) {
				double vec =
					sqrt(pow(V[i].x - V[j].x, 2) + pow(V[i].y - V[j].y, 2));
				if (vec > vec_max) {
					vec_max = vec;
					index_i = i;
					index_j = j;
				}
			}
		}
		return std::make_pair(vec_max, std::make_pair(index_i, index_j));
	}

	double polygonArea(const Vertices& V) {
		/*
		An implementation of the shoelace algorithm, first described by Albrecht
		Ludwig Friedrich Meister, which is used to calculate the area of a
		polygon.
		*/

		const unsigned long N = V.size();
		double out = 0.;
		for (unsigned long n = 0; n < N; n++) {
			out += V[n].x * V[(n + 1) % N].y - V[n].y * V[(n + 1) % N].x;
		}
		return abs(out) / 2.;
	}

}}