/*
Functions for producing group theoretic transformations.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <math.h>
#include <utility>

// src
#include "../types.hpp"
#include "./polygon_properties.hpp"
using namespace kac_core::types;

namespace kac_core { namespace geometry {

	Vertices convexNormalisation(Vertices V) {
		/*
		This algorithm produces an identity polygon for each unique polygon
		given as input. This method normalises an input polygon to the unit
		interval such that x ∈ [0, 1] && y ∈ [0, 1], reducing each input polygon
		by isometric and similarity transformations. This is achieved by first
		enforcing that the vertices of a polygon are ordered clockwise. Then,
		the largest vector is used to determine the lower and upper bounds
		across the x-axis. Next, the polygon is split into quadrants, the
		largest of whose area determines the rotation/reflection of the polygon.
		Finally, the points are normalised, and ordered such that
		V[0] = [0., y].
		*/

		// enforce that each polygon is clockwise
		// reverse the polygon if the vertices are anti-clockwise
		if ((V[1].x - V[0].x) * (V[2].y - V[1].y)
				- (V[2].x - V[1].x) * (V[1].y - V[0].y)
			> 0) {
			std::reverse(V.begin(), V.end());
		}

		// orient largest vector across x-axis
		// determine largest vector
		std::pair<double, std::pair<int, int>> LV = largestVector(V);
		// shift midpoint of the largest vector to origin
		double x_shift = (V[LV.second.first].x + V[LV.second.second].x) / 2;
		double y_shift = (V[LV.second.first].y + V[LV.second.second].y) / 2;
		for (unsigned long n = 0; n < V.size(); n++) {
			V[n].x -= x_shift;
			V[n].y -= y_shift;
		}
		// rotate around midpoint such that largest_vec is vertical
		double theta = V[LV.second.first].theta();
		double cos_theta = cos(theta);
		double sin_theta = sin(theta);
		for (unsigned long n = 0; n < V.size(); n++) {
			V[n] = Point(
				V[n].x * cos_theta + V[n].y * sin_theta,
				-V[n].x * sin_theta + V[n].y * cos_theta
			);
		}

		// normalise to unit interval
		// first find minmax in both x & y
		Matrix_1D X;
		Matrix_1D Y;
		for (unsigned long n = 0; n < V.size(); n++) {
			X.push_back(V[n].x);
			Y.push_back(V[n].y);
		}
		auto x_min_max = std::minmax_element(begin(X), end(X));
		auto y_min_max = std::minmax_element(begin(Y), end(Y));
		// center along y-axis (reuse variable y_shift)
		y_shift = (*y_min_max.first + *y_min_max.second) / 2;
		for (unsigned long n = 0; n < V.size(); n++) { V[n].y -= y_shift; }
		// find overall v_min
		*y_min_max.first -= y_shift;
		*y_min_max.second -= y_shift;
		double v_min = *x_min_max.first < *y_min_max.first ? *x_min_max.first
														   : *y_min_max.first;
		double v_d = (*x_min_max.second > *y_min_max.second ? *x_min_max.second
															: *y_min_max.second)
			- v_min;
		// normalise
		for (unsigned long n = 0; n < V.size(); n++) {
			V[n].x = (V[n].x - v_min) / v_d;
			V[n].y = (V[n].y - v_min) / v_d;
		}

		// position x = 0. at V[0]
		unsigned long n_shift = 0;
		for (unsigned long n = 0; n < V.size(); n++) {
			if (V[n].x == 0.) {
				n_shift = n;
				break;
			}
		}
		std::rotate(V.begin(), V.begin() + n_shift, V.end());

		return V;
	}

	// Point squareToCircle(const Point& p) {
	// 	/*
	// 	Project a point from the Euclidean plane onto a circular
	// lattice. This 	projection is neither equiareal nor conformal,
	// but is a compromise 	between the two.

	// 	See here:
	// 		Fong, C. (2014). Analytical methods for squaring the disc.
	// 		http://arxiv.org/abs/1509.06344
	// 	*/

	// 	return Point(
	// 		p.x * sqrt(1 - p.y * p.y / 2), p.y * sqrt(1 - p.x * p.x / 2)
	// 	);
	// }

}}