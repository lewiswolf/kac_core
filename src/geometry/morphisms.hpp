/*
Functions for producing group theoretic transformations.
*/

#pragma once

// core
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
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
		double x_min = P[0].x;
		double x_max = P[0].x;
		double y_min = P[0].y;
		double y_max = P[0].y;
		for (T::Point& p: P) {
			x_min = std::min(x_min, p.x);
			x_max = std::max(x_max, p.x);
			y_min = std::min(y_min, p.y);
			y_max = std::max(y_max, p.y);
		}
		// center along x and y axes
		const double x_shift = (x_min + x_max) * 0.5;
		const double y_shift = (y_min + y_max) * 0.5;
		for (T::Point& p: P) {
			p.x -= x_shift;
			p.y -= y_shift;
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
			for (T::Point& p: P) {
				p.x = 2 * (p.x - v_min) / v_d - 1;
				p.y = 2 * (p.y - v_min) / v_d - 1;
			}
		} else {
			for (T::Point& p: P) {
				p.x = (p.x - v_min) / v_d;
				p.y = (p.y - v_min) / v_d;
			}
		}
		// enforce that each polygon is anti-clockwise
		// reverse the polygon if the vertices are clockwise
		if (polygonArea(P) < 0.) {
			std::reverse(P.begin(), P.end());
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
		if (polygonArea(P) > 0.) {
			std::reverse(P.begin(), P.end());
		}
		// orient largest vector across x-axis
		const std::size_t N = P.size();
		// determine largest vector
		const auto& [_, indices] = largestVector(P);
		// shift midpoint of the largest vector to origin
		const double x_shift = (P[indices.first].x + P[indices.second].x) * 0.5;
		const double y_shift = (P[indices.first].y + P[indices.second].y) * 0.5;
		for (T::Point& p: P) {
			p.x -= x_shift;
			p.y -= y_shift;
		}
		// rotate around midpoint such that largest_vec is horizontal
		const double theta = P[indices.first].theta();
		const double cos_theta = std::cos(theta);
		const double sin_theta = std::sin(theta);
		for (T::Point& p: P) {
			p = T::Point(p.x * cos_theta + p.y * sin_theta, -p.x * sin_theta + p.y * cos_theta);
		}
		// find area of each cartesian quadrant and position the largest in the positive x and y to
		// remove isometric transformations.
		std::array<double, 4> quadAreas = {0., 0., 0., 0.};
		// determine which quadrant each point is in
		const auto whichQuad = [](const T::Point& p) {
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
		const auto triangleArea = [](const T::Point& a, const T::Point& b) {
			return abs(b.y * a.x - b.x * a.y) * 0.5;
		};
		// loop over points and sum quadrant areas
		for (std::size_t n = 0; n < N; n++) {
			const T::Point a = P[n];
			const T::Point b = P[(n + 1) % N];
			const short quad_a = whichQuad(a);
			const short quad_b = whichQuad(b);
			if (quad_a == quad_b) {
				// if the points lie in the same quadrant, add the triangular area to that quadrant.
				quadAreas[quad_a] += triangleArea(a, b);
			} else if (((quad_b + 4) - quad_a) % 4 == 1) {
				// if the points are in neighbouring quadrants, find the axis intersection between
				// them and update both quadrants.
				T::Point c = T::Point(0., 0.);
				if ((a.x * b.x) < 0.) {
					// crosses y axis
					c.x = 0.;
					c.y = a.y - (b.y - a.y) / (b.x - a.x) * a.x;
				} else {
					// crosses x axis
					c.x = a.x - (b.x - a.x) / (b.y - a.y) * a.y;
					c.y = 0.;
				}
				quadAreas[quad_a] += triangleArea(a, c);
				quadAreas[quad_b] += triangleArea(c, b);
			} else if (((quad_b + 4) - quad_a) % 4 == 2) {
				// if the points lie across three quadrants, update the two quadrants containing the
				// points in question.
				T::Point c = T::Point(0, a.y - (b.y - a.y) / (b.x - a.x) * a.x);
				T::Point d = T::Point(a.x - (b.x - a.x) / (b.y - a.y) * a.y, 0);
				if (std::hypot(a.x - c.x, a.y - c.y) < std::hypot(a.x - d.x, a.y - d.y)) {
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
				for (std::size_t n = 0; n < N; n++) { P[n].y *= -1.; }
				std::reverse(P.begin(), P.end());
				break;
			case 2:
				for (std::size_t n = 0; n < N; n++) {
					P[n] = T::Point(P[n].x *= -1., P[n].y *= -1.);
				}
				break;
			case 3:
				for (std::size_t n = 0; n < N; n++) { P[n].x *= -1.; }
				std::reverse(P.begin(), P.end());
				break;
		}
		// normalise
		P = normalisePolygon(P, signed_norm);
		// position x = -1 : 0 at P[0]
		std::size_t n_shift = 0;
		if (signed_norm) {
			for (std::size_t n = 0; n < N; n++) {
				if (P[n].x == -1.) {
					n_shift = n;
					break;
				}
			}
		} else {
			for (std::size_t n = 0; n < N; n++) {
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
		if (polygonArea(P) > 0.) {
			std::reverse(P.begin(), P.end());
		}
		// orient largest vector across x-axis
		const std::size_t N = P.size();
		// determine largest vector
		const auto& [_, indices] = largestVector(P);
		// shift midpoint of the largest vector to origin
		const double x_shift = (P[indices.first].x + P[indices.second].x) * 0.5;
		const double y_shift = (P[indices.first].y + P[indices.second].y) * 0.5;
		for (T::Point& p: P) {
			p.x -= x_shift;
			p.y -= y_shift;
		}
		// rotate around midpoint such that largest_vec is horizontal
		const double theta = P[indices.first].theta();
		const double cos_theta = std::cos(theta);
		const double sin_theta = std::sin(theta);
		for (T::Point& p: P) {
			p = T::Point(p.x * cos_theta + p.y * sin_theta, -p.x * sin_theta + p.y * cos_theta);
		}
		// normalise
		P = normalisePolygon(P, signed_norm);
		// position x = -1 : 0 at P[0]
		std::size_t n_shift = 0;
		if (signed_norm) {
			for (std::size_t n = 0; n < N; n++) {
				if (P[n].x == -1.) {
					n_shift = n;
					break;
				}
			}
		} else {
			for (std::size_t n = 0; n < N; n++) {
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

		const double scale = std::sqrt(std::abs(a) / std::abs(polygonArea(P)));
		const T::Point centroid = polygonCentroid(P);
		for (T::Point& p: P) {
			p.x = centroid.x + scale * (p.x - centroid.x);
			p.y = centroid.y + scale * (p.y - centroid.y);
		}
		// enforce signed area
		if ((polygonArea(P) > 0.) != (a > 0.)) {
			std::reverse(P.begin(), P.end());
		}
		return P;
	}

}
