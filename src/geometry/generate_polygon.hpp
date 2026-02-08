/*
Functions for generating polygons.
*/

#pragma once

// core
#include <algorithm>	// generate, min, max, shuffle, sort
#include <cmath>		// abs
#include <limits>		// numeric_limits
#include <math.h>		// abs
#include <numbers>		// pi
#include <random>		// default_random_engine
#include <time.h>		// time

// src
#include "../types.hpp"
#include "./lines.hpp"
namespace T = kac_core::types;

// random number generators
static std::default_random_engine random_engine(time(0l));
static std::uniform_real_distribution<double> bipolar_distribution(-1., 1.);
static std::uniform_int_distribution<long>
	uniform_sequence(std::numeric_limits<long>::min(), std::numeric_limits<long>::max());
// numerical constants
constexpr double infinity = std::numeric_limits<double>::infinity();

namespace kac_core::geometry {

	inline T::Polygon generateConvexPolygon(const unsigned long& N, const time_t& seed = 0l) {
		/*
		Generate convex shapes according to Pavel Valtr's 1995 algorithm.
		Adapted from Sander Verdonschot's Java version, found here:
		https://cglab.ca/~sander/misc/ConvexGeneration/ValtrAlgorithm.java
		input:
			N = the number of vertices
			seed? = the seed for the random number generators
		output:
			P = a convex polygon of N random vertices
		*/

		// initialise variables
		T::Polygon P;
		std::vector<double> X(N);
		std::vector<double> Y(N);
		std::vector<double> X_rand(N);
		std::vector<double> Y_rand(N);
		unsigned long last_true = 0;
		unsigned long last_false = 0;
		// initialise and sort random coordinates
		if (seed != 0l) {
			random_engine.seed(seed);
		}
		std::generate(X_rand.begin(), X_rand.end(), []() {
			return bipolar_distribution(random_engine);
		});
		std::generate(Y_rand.begin(), Y_rand.end(), []() {
			return bipolar_distribution(random_engine);
		});
		std::sort(X_rand.begin(), X_rand.end());
		std::sort(Y_rand.begin(), Y_rand.end());
		// divide the interior points into two chains
		for (unsigned long n = 1; n < N; n++) {
			if (n != N - 1) {
				if (uniform_sequence(random_engine) % 2 == 1) {
					X[n] = X_rand[n] - X_rand[last_true];
					Y[n] = Y_rand[n] - Y_rand[last_true];
					last_true = n;
				} else {
					X[n] = X_rand[last_false] - X_rand[n];
					Y[n] = Y_rand[last_false] - Y_rand[n];
					last_false = n;
				}
			} else {
				X[0] = X_rand[n] - X_rand[last_true];
				Y[0] = Y_rand[n] - Y_rand[last_true];
				X[n] = X_rand[last_false] - X_rand[n];
				Y[n] = Y_rand[last_false] - Y_rand[n];
			}
		}
		// randomly combine x and y
		shuffle(Y.begin(), Y.end(), random_engine);
		for (unsigned long n = 0; n < N; n++) { P.push_back(T::Point(X[n], Y[n])); }
		// sort by polar angle
		sort(P.begin(), P.end(), [](T::Point& p1, T::Point& p2) {
			return p1.theta() < p2.theta();
		});
		// arrange points end to end to form a polygon
		double x_min, x_max, y_min, y_max = 0;
		double x = 0.0;
		double y = 0.0;
		for (unsigned long n = 0; n < N; n++) {
			T::Point p = T::Point(x, y);
			x += P[n].x;
			y += P[n].y;
			P[n] = p;
			x_min = std::min(P[n].x, x_min);
			x_max = std::max(P[n].x, x_max);
			y_min = std::min(P[n].y, y_min);
			y_max = std::max(P[n].y, y_max);
		}
		// center around origin
		double x_shift = ((x_max - x_min) * 0.5) - x_max;
		double y_shift = ((y_max - y_min) * 0.5) - y_max;
		for (unsigned long n = 0; n < N; n++) {
			P[n].x += x_shift;
			P[n].y += y_shift;
		}
		return P;
	}

	inline T::Polygon generateIrregularStar(const std::size_t& N, const time_t& seed = 0l) {
		/*
		This is a fast method for generating concave polygons, particularly with a large number of
		vertices. This approach generates polygons by ordering a series of random points around a
		centre point. As a result, not all possible simple polygons are generated this way.
		input:
			N = the number of vertices
			seed? = the seed for the random number generators
		output:
			P = an irregular star of N random vertices
		*/

		// initialise variables
		T::Matrix_1D X(N, 0.);
		T::Matrix_1D Y(N, 0.);
		T::Polygon P(N, T::Point(0., 0.));
		// initialise and sort random coordinates
		if (seed != 0l) {
			random_engine.seed(seed);
		}
		// first find minmax in both x & y
		for (std::size_t n = 0; n < N; n++) {
			X[n] = bipolar_distribution(random_engine);
			Y[n] = bipolar_distribution(random_engine);
		}
		// center along x and y axes
		double x_min = infinity;
		double x_max = -infinity;
		double y_min = infinity;
		double y_max = -infinity;
		for (std::size_t n = 0; n < N; n++) {
			x_min = std::min(X[n], x_min);
			x_max = std::max(X[n], x_max);
			y_min = std::min(Y[n], y_min);
			y_max = std::max(Y[n], y_max);
		}
		const double x_shift = (x_min + x_max) * 0.5;
		const double y_shift = (y_min + y_max) * 0.5;
		for (std::size_t n = 0; n < N; n++) {
			P[n].x = X[n] - x_shift;
			P[n].y = Y[n] - y_shift;
			// normalise to unit circle
			P[n].x = P[n].r() * std::cos(P[n].theta()) / std::numbers::sqrt2;
			P[n].y = P[n].r() * std::sin(P[n].theta()) / std::numbers::sqrt2;
		}
		// sort by polar angle
		sort(P.begin(), P.end(), [](T::Point& a, T::Point& b) { return a.theta() < b.theta(); });
		return P;
	}

	inline T::Polygon generatePolygon(const std::size_t& N, const time_t& seed = 0l) {
		/*
		This algorithm is based on a method of eliminating self-intersections in a polygon by
		using the Lin and Kerningham '2-opt' moves. Such a move eliminates an intersection between
		two edges by reversing the order of the vertices between the edges. Intersecting edges are
		detected using a simple sweep through the vertices and then one intersection is chosen at
		random to eliminate after each sweep.
		https://doc.cgal.org/latest/Generator/group__PkgGeneratorsRef.html#gaa8cb58e4cc9ab9e225808799b1a61174
		van Leeuwen, J., & Schoone, A. A. (1982). Untangling a traveling salesman tour in the plane.
		input:
			N = the number of vertices
			seed? = the seed for the random number generators
		output:
			P = a concave polygon of N random vertices
		*/

		// initialise variables
		T::Polygon P(N, T::Point(0., 0.));
		// initialise and sort random coordinates
		if (seed != 0l) {
			random_engine.seed(seed);
		}
		for (std::size_t n = 0; n < N; n++) {
			P[n].x = bipolar_distribution(random_engine);
			P[n].y = bipolar_distribution(random_engine);
		}
		// 2 opt loop
		std::vector<std::pair<std::size_t, std::size_t>> indices;
		bool intersections = true;
		while (intersections) {
			bool restarted = false;
			for (std::size_t i = 0; i < N - 1 && !restarted; i++) {
				for (std::size_t j = i + 1; j < N && !restarted; j++) {
					// collect indices of lines which should be crossed
					auto intersection_type =
						lineIntersection(T::Line(P[i], P[i + 1]), T::Line(P[j], P[(j + 1) % N]))
							.first;
					if (intersection_type == IntersectionType::None
						|| intersection_type == IntersectionType::Vertex) {
						continue;
					} else if (intersection_type == IntersectionType::Intersect) {
						indices.emplace_back(i + 1, j + 1);
					} else if (intersection_type == IntersectionType::Branch) {
						indices.emplace_back(i, j);
					} else if (intersection_type == IntersectionType::Colinear) {
						std::reverse(
							P.begin() + i + (P[i].x < P[i + 1].x ? 0 : 1),
							P.begin() + j + (P[j].x > P[j + 1].x ? 0 : 1)
						);
						// restart loop
						indices.clear();
						restarted = true;
					}
				}
			}
			if (restarted) {
				// goto replacement
				continue;
			}
			if (indices.empty()) {
				// close loop
				intersections = false;
			} else {
				std::pair<std::size_t, std::size_t> swap =
					indices[abs(uniform_sequence(random_engine)) % indices.size()];
				std::reverse(P.begin() + swap.first, P.begin() + swap.second);
				// restart loop
				indices.clear();
			}
		}
		return P;
	}

	inline T::Polygon generateRegularPolygon(const std::size_t& N) {
		/*
		Generate a N-sided regular polygon.
		*/

		T::Polygon P(N, T::Point(0., 0.));
		double theta = 0.;
		const double d_theta = 2. * std::numbers::pi / N;
		for (std::size_t n = 0; n < N; n++) {
			P[n].x = std::cos(theta);
			P[n].y = std::sin(theta);
			theta += d_theta;
		}
		return P;
	}

	inline T::Polygon generateUnitRectangle(const double& epsilon) {
		/*
		Define a rectangle with unit area and an aspect ration epsilon.
		*/

		double x = 0.;
		double y = infinity;
		if (epsilon != 0.) {
			x = 0.5 * epsilon;
			y = 0.5 / epsilon;
		}
		return T::Polygon({T::Point(x, y), T::Point(-x, y), T::Point(-x, -y), T::Point(x, -y)});
	}

	inline T::Polygon generateUnitTriangle(const double& r, const double& theta) {
		/*
		Define a triangle with unit area using a polar coordinate point mapped onto a lens.
		See Guy, R. K. (1993). There are three times as many obtuse-angled triangles as there are
		acute-angled ones.
		*/

		const double cos_theta = std::cos(theta);
		const double sin_theta = std::sin(theta);
		// map unit disk -> radial boundary of the lens
		const double rho =
			r * (std::sqrt(1. - 0.25 * sin_theta * sin_theta) - 0.5 * std::abs(cos_theta));
		// calculate the mapped point
		T::Point p = T::Point(rho * cos_theta, rho * sin_theta);
		// enforce unit area
		const double height = std::abs(p.y);
		if (height == 0) {
			return T::Polygon({T::Point(-infinity, 0.), p, T::Point(infinity, 0.)});
		} else {
			const double scale = 1. / std::sqrt(0.5 * height);
			p.x *= scale;
			p.y *= scale;
			if (p.y > 0.) {
				return T::Polygon({T::Point(-0.5 * scale, 0.), T::Point(0.5 * scale, 0.), p});
			} else {
				return T::Polygon({T::Point(0.5 * scale, 0.), T::Point(-0.5 * scale, 0.), p});
			}
		}
	}

}
