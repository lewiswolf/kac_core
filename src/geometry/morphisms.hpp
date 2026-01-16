/*
Functions for producing group theoretic transformations.
*/

#pragma once

// core
#include <algorithm>
#include <array>
#include <iterator>
#include <math.h>
#include <utility>

// src
#include "../types.hpp"
#include "./polygon_properties.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline T::Polygon normalisePolygon(T::Polygon P, const bool& signed_norm = false) {
		/*
		This function takes a polygon, centers it across the x and y axis, then normalises the
		vertices to the unit interval ℝ^2.
		*/

		// first find minmax in both x & y
		const std::size_t N = P.size();
		double x_min = P[0].x, x_max = P[0].x;
		double y_min = P[0].y, y_max = P[0].y;
		for (std::size_t n = 1; n < N; ++n) {
			x_min = std::min(x_min, P[n].x);
			x_max = std::max(x_max, P[n].x);
			y_min = std::min(y_min, P[n].y);
			y_max = std::max(y_max, P[n].y);
		}
		// center along x and y axes
		const double x_shift = (x_min + x_max) / 2;
		const double y_shift = (y_min + y_max) / 2;
		for (std::size_t n = 0; n < N; n++) {
			P[n].x -= x_shift;
			P[n].y -= y_shift;
		}
		x_min -= x_shift;
		x_max -= x_shift;
		y_min -= y_shift;
		y_max -= y_shift;
		// find v_min and v_d (v_d = v_max - v_min)
		const double v_min = std::min(x_min, y_min);
		const double v_d = std::max(x_max, y_max) - v_min;
		// normalise
		if (signed_norm) {
			for (std::size_t n = 0; n < N; n++) {
				P[n].x = 2 * (P[n].x - v_min) / v_d - 1;
				P[n].y = 2 * (P[n].y - v_min) / v_d - 1;
			}
		} else {
			for (std::size_t n = 0; n < N; n++) {
				P[n].x = (P[n].x - v_min) / v_d;
				P[n].y = (P[n].y - v_min) / v_d;
			}
		}
		return P;
	}

	inline T::Polygon normaliseConvexPolygon(T::Polygon P, const bool& signed_norm = false) {
		/*
		This algorithm produces an identity polygon for each unique polygon given as input. This
		method normalises an input polygon to the unit interval such that x ∈ [0, 1] && y ∈ [0, 1],
		reducing each input polygon by isometric and similarity transformations. This is achieved by
		first enforcing that the vertices of a polygon are ordered clockwise. Then, the largest
		vector is used to determine the lower and upper bounds across the x-axis. Next, the polygon
		is split into quadrants, the largest of whose area determines the rotation/reflection of the
		polygon. Finally, the points are normalised, and ordered such that P[0] = [0., y].
		*/

		// enforce that each polygon is clockwise
		// reverse the polygon if the vertices are anti-clockwise
		if ((P[1].x - P[0].x) * (P[2].y - P[1].y) - (P[2].x - P[1].x) * (P[1].y - P[0].y) > 0) {
			std::reverse(P.begin(), P.end());
		}
		// orient largest vector across x-axis
		// determine largest vector
		std::pair<double, std::pair<long, long>> LV = largestVector(P);
		// shift midpoint of the largest vector to origin
		double x_shift = (P[LV.second.first].x + P[LV.second.second].x) / 2;
		double y_shift = (P[LV.second.first].y + P[LV.second.second].y) / 2;
		const unsigned long N = P.size();
		for (unsigned long n = 0; n < N; n++) {
			P[n].x -= x_shift;
			P[n].y -= y_shift;
		}
		// rotate around midpoint such that largest_vec is horizontal
		double theta = P[LV.second.first].theta();
		double cos_theta = cos(theta);
		double sin_theta = sin(theta);
		for (unsigned long n = 0; n < N; n++) {
			P[n] = T::Point(
				P[n].x * cos_theta + P[n].y * sin_theta, -P[n].x * sin_theta + P[n].y * cos_theta
			);
		}
		// find area of each cartesian quadrant and position the largest in the positive x and y to
		// remove isometric transformations.
		std::array<double, 4> quadAreas = {0., 0., 0., 0.};
		// determine which quadrant each point is in
		auto whichQuad = [](const T::Point& p) {
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
		auto triangleArea = [](const T::Point& a, const T::Point& b) {
			return abs(b.y * a.x - b.x * a.y) / 2;
		};
		// loop over points and sum quadrant areas
		for (unsigned long n = 0; n < N; n++) {
			T::Point a = P[n];
			T::Point b = P[(n + 1) % N];
			short quad_a = whichQuad(a);
			short quad_b = whichQuad(b);
			if (quad_a == quad_b) {
				// if the points lie in the same quadrant, add the triangular area to that quadrant.
				quadAreas[quad_a] += triangleArea(a, b);
			} else if (((quad_b + 4) - quad_a) % 4 == 1) {
				// if the points are in neighbouring quadrants, find the axis intersection between
				// them and update both quadrants.
				T::Point c;
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
				// if the points lie across three quadrants, update the two quadrants containing the
				// points in question.
				T::Point c = T::Point(0, a.y - (b.y - a.y) / (b.x - a.x) * a.x);
				T::Point d = T::Point(a.x - (b.x - a.x) / (b.y - a.y) * a.y, 0);
				if (hypot(a.x - c.x, a.y - c.y) < hypot(a.x - d.x, a.y - d.y)) {
					quadAreas[quad_a] += triangleArea(a, c);
					// quadAreas[(quad_a + 1) % 4] += triangleArea(c, d);
					quadAreas[quad_b] += triangleArea(d, b);
				} else {
					quadAreas[quad_a] += triangleArea(a, d);
					// quadAreas[(quad_a + 1) % 4] += triangleArea(d, c);
					quadAreas[quad_b] += triangleArea(c, b);
				}
			}
		}
		// reflect initial polygon given the largest quadrant
		switch (static_cast<short>(
			std::distance(quadAreas.begin(), std::max_element(quadAreas.begin(), quadAreas.end()))
		)) {
			case 0:
				break;
			case 1:
				for (unsigned long n = 0; n < N; n++) { P[n].y *= -1.; }
				std::reverse(P.begin(), P.end());
				break;
			case 2:
				for (unsigned long n = 0; n < N; n++) {
					P[n] = T::Point(P[n].x *= -1., P[n].y *= -1.);
				}
				break;
			case 3:
				for (unsigned long n = 0; n < N; n++) { P[n].x *= -1.; }
				std::reverse(P.begin(), P.end());
				break;
		}
		// normalise
		P = normalisePolygon(P, signed_norm);
		// position x = -1 : 0 at P[0]
		unsigned long n_shift = 0;
		if (signed_norm) {
			for (unsigned long n = 0; n < N; n++) {
				if (P[n].x == -1.) {
					n_shift = n;
					break;
				}
			}
		} else {
			for (unsigned long n = 0; n < N; n++) {
				if (P[n].x == -0.) {
					n_shift = n;
					break;
				}
			}
		}
		std::rotate(P.begin(), P.begin() + n_shift, P.end());
		return P;
	}

	inline T::Polygon normaliseSimplePolygon(T::Polygon P, const bool& signed_norm = false) {
		/*
		This algorithm performs general normalisation rotations to ensure uniqueness, however it is
		not comprehensive for all simple geometric transformations.
		*/

		// enforce that each polygon is clockwise
		// reverse the polygon if the vertices are anti-clockwise
		if ((P[1].x - P[0].x) * (P[2].y - P[1].y) - (P[2].x - P[1].x) * (P[1].y - P[0].y) > 0) {
			std::reverse(P.begin(), P.end());
		}
		// orient largest vector across x-axis
		// determine largest vector
		std::pair<double, std::pair<long, long>> LV = largestVector(P);
		// shift midpoint of the largest vector to origin
		double x_shift = (P[LV.second.first].x + P[LV.second.second].x) / 2;
		double y_shift = (P[LV.second.first].y + P[LV.second.second].y) / 2;
		const unsigned long N = P.size();
		for (unsigned long n = 0; n < N; n++) {
			P[n].x -= x_shift;
			P[n].y -= y_shift;
		}
		// rotate around midpoint such that largest_vec is horizontal
		double theta = P[LV.second.first].theta();
		double cos_theta = cos(theta);
		double sin_theta = sin(theta);
		for (unsigned long n = 0; n < N; n++) {
			P[n] = T::Point(
				P[n].x * cos_theta + P[n].y * sin_theta, -P[n].x * sin_theta + P[n].y * cos_theta
			);
		}
		// normalise
		P = normalisePolygon(P, signed_norm);
		// position x = -1 : 0 at P[0]
		unsigned long n_shift = 0;
		if (signed_norm) {
			for (unsigned long n = 0; n < N; n++) {
				if (P[n].x == -1.) {
					n_shift = n;
					break;
				}
			}
		} else {
			for (unsigned long n = 0; n < N; n++) {
				if (P[n].x == -0.) {
					n_shift = n;
					break;
				}
			}
		}
		std::rotate(P.begin(), P.begin() + n_shift, P.end());
		return P;
	}

	inline T::Polygon scalePolygonByArea(T::Polygon P, const double& a) {
		/*
		Scale a polygon by area, whilst preserving angle and distance relationships between
		vertices.
		*/

		double scale = std::sqrt(std::abs(a) / polygonArea(P));
		T::Point centroid = polygonCentroid(P);
		for (std::size_t n = 0; n < P.size(); n++) {
			P[n].x = centroid.x + scale * (P[n].x - centroid.x);
			P[n].y = centroid.y + scale * (P[n].y - centroid.y);
		}
		return P;
	}

}
