/*
Functions for producing group theoretic transformations.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <iterator>
#include <math.h>
#include <utility>

// src
#include "../types.hpp"
#include "./polygon_properties.hpp"
using namespace kac_core::types;

namespace kac_core { namespace geometry {

	Vertices normalisePolygon(Vertices V) {
		/*
		This function takes a polygon, centers it across the x and y axis, then
		normalises the vertices to the unit interval ℝ^2.
		*/

		// first find minmax in both x & y
		Matrix_1D X;
		Matrix_1D Y;
		for (unsigned long n = 0; n < V.size(); n++) {
			X.push_back(V[n].x);
			Y.push_back(V[n].y);
		}
		auto x_min_max = std::minmax_element(begin(X), end(X));
		auto y_min_max = std::minmax_element(begin(Y), end(Y));
		// center along x and y axes
		double x_shift = (*x_min_max.first + *x_min_max.second) / 2;
		double y_shift = (*y_min_max.first + *y_min_max.second) / 2;
		for (unsigned long n = 0; n < V.size(); n++) {
			V[n].x -= x_shift;
			V[n].y -= y_shift;
		}
		*x_min_max.first -= x_shift;
		*x_min_max.second -= x_shift;
		*y_min_max.first -= y_shift;
		*y_min_max.second -= y_shift;
		// find v_min and v_d (v_d = v_max - v_min)
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
		return V;
	}

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
		// find area of each cartesian quadrant and position the largest in
		// the positive x and y to remove isometric transformations.
		std::array<double, 4> quadAreas = {0., 0., 0., 0.};
		// determine which quadrant each point is in
		auto whichQuad = [](const Point& p) {
			if (p.x >= 0. && p.y > 0.) {
				return 0;
			} else if (p.x > 0. && p.y <= 0.) {
				return 1;
			} else if (p.x <= 0. && p.y < 0.) {
				return 2;
			} else {
				return 3;
			}
		};
		// area of a triangle, simplified due to point c = [0., 0.]
		auto triangleArea = [](const Point& a, const Point& b) {
			return fabs(b.y * a.x - b.x * a.y) / 2;
		};
		// loop over points and sum quadrant areas
		for (unsigned long n = 0; n < V.size(); n++) {
			Point a = V[n];
			Point b = V[(n + 1) % V.size()];
			int quad_a = whichQuad(a);
			int quad_b = whichQuad(b);
			if (quad_a == quad_b) {
				// if the points lie in the same quadrant, add the triangular
				// area to that quadrant.
				quadAreas[quad_a] += triangleArea(a, b);
			} else if (((quad_b + 4) - quad_a) % 4 == 1) {
				// if the points are in neighbouring quadrants, find the axis
				// intersection between them and update both quadrants.
				Point c;
				if ((a.x * b.x) < 0) {
					// crosses y axis
					c.x = 0;
					c.y = a.y - (b.y - a.y) / (b.x - a.x) * a.x;
				} else {
					// crosses x axis
					c.x = a.x - (b.x - a.x) / (b.y - a.y) * a.y;
					c.y = 0;
				}
				quadAreas[quad_a] += triangleArea(a, c);
				quadAreas[quad_b] += triangleArea(c, b);
			} else if (((quad_b + 4) - quad_a) % 4 == 2) {
				// if the points lie across three quadrants, update the two
				// quadrants as well as the one in between.
				Point c = Point(0, a.y - (b.y - a.y) / (b.x - a.x) * a.x);
				Point d = Point(a.x - (b.x - a.x) / (b.y - a.y) * a.y, 0);
				if (sqrt(pow(a.x - c.x, 2) + pow(a.y - c.y, 2))
					< sqrt(pow(a.x - d.x, 2) + pow(a.y - d.y, 2))) {
					quadAreas[quad_a] += triangleArea(a, c);
					quadAreas[quad_a + 1 % 4] += triangleArea(c, d);
					quadAreas[quad_b] += triangleArea(d, b);
				} else {
					quadAreas[quad_a] += triangleArea(a, d);
					quadAreas[quad_a + 1 % 4] += triangleArea(d, c);
					quadAreas[quad_b] += triangleArea(c, b);
				}
			}
		}
		// reflect initial polygon given the largest quadrant
		switch (static_cast<int>(std::distance(
			quadAreas.begin(),
			std::max_element(quadAreas.begin(), quadAreas.end())
		))) {
			case 0:
				break;
			case 1:
				for (unsigned long n = 0; n < V.size(); n++) { V[n].y *= -1.; }
				std::reverse(V.begin(), V.end());
				break;
			case 2:
				for (unsigned long n = 0; n < V.size(); n++) {
					V[n] = Point(V[n].x *= -1., V[n].y *= -1.);
				}
				break;
			case 3:
				for (unsigned long n = 0; n < V.size(); n++) { V[n].x *= -1.; }
				std::reverse(V.begin(), V.end());
				break;
		}
		// normalise
		V = normalisePolygon(V);
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

}}