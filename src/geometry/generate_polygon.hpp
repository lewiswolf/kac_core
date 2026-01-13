/*
Functions for generating polygons.
*/

#pragma once

// core
#include <algorithm>	// generate, min, max, shuffle, sort
#include <limits>		// numeric_limits
#include <math.h>		// abs
#include <random>		// default_random_engine
#include <string>		// string
#include <time.h>		// time

// src
#include "../types.hpp"
#include "./lines.hpp"
namespace T = kac_core::types;

// random number generators
static std::default_random_engine random_engine(time(0l));
static std::uniform_real_distribution<double> uniform_distribution(-1., 1.);
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
			return uniform_distribution(random_engine);
		});
		std::generate(Y_rand.begin(), Y_rand.end(), []() {
			return uniform_distribution(random_engine);
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
		double x_shift = ((x_max - x_min) / 2.0) - x_max;
		double y_shift = ((y_max - y_min) / 2.0) - y_max;
		for (unsigned long n = 0; n < N; n++) {
			P[n].x += x_shift;
			P[n].y += y_shift;
		}
		return P;
	}

	inline T::Polygon generateIrregularStar(const unsigned long& N, const time_t& seed = 0l) {
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
		T::Matrix_1D X;
		T::Matrix_1D Y;
		T::Polygon P;
		// initialise and sort random coordinates
		if (seed != 0l) {
			random_engine.seed(seed);
		}
		// first find minmax in both x & y
		for (unsigned long n = 0; n < N; n++) {
			X.push_back(uniform_distribution(random_engine));
			Y.push_back(uniform_distribution(random_engine));
		}
		auto x_min_max = std::minmax_element(begin(X), end(X));
		auto y_min_max = std::minmax_element(begin(Y), end(Y));
		// center along x and y axes
		double x_shift = (*x_min_max.first + *x_min_max.second) / 2;
		double y_shift = (*y_min_max.first + *y_min_max.second) / 2;
		for (unsigned long n = 0; n < N; n++) {
			P.push_back(T::Point(X[n] -= x_shift, Y[n] -= y_shift));
		}
		// sort by polar angle
		sort(P.begin(), P.end(), [](T::Point& a, T::Point& b) { return a.theta() < b.theta(); });
		return P;
	}

	inline T::Polygon generatePolygon(const unsigned long& N, const time_t& seed = 0l) {
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
		T::Polygon P;
		// initialise and sort random coordinates
		if (seed != 0l) {
			random_engine.seed(seed);
		}
		for (unsigned long n = 0; n < N; n++) {
			P.push_back(
				T::Point(uniform_distribution(random_engine), uniform_distribution(random_engine))
			);
		}
		// 2 opt loop
		std::vector<std::pair<long, long>> indices;
		std::string intersection_type = "";
		bool intersections = true;
		while (intersections) {
			bool restarted = false;
			for (unsigned long i = 0; i < N - 2 && !restarted; i++) {
				for (unsigned long j = i + 1; j < N && !restarted; j++) {
					// collect indices of lines which should be crossed
					intersection_type =
						lineIntersection(T::Line(P[i], P[i + 1]), T::Line(P[j], P[(j + 1) % N]))
							.first;
					if (intersection_type == "none" || intersection_type == "vertex") {
						continue;
					} else if (intersection_type == "intersect") {
						indices.emplace_back(i + 1, j + 1);
					} else if (intersection_type == "adjacent") {
						indices.emplace_back(i, j);
					} else if (intersection_type == "colinear") {
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
				std::pair<long, long> swap =
					indices[abs(uniform_sequence(random_engine)) % indices.size()];
				std::reverse(P.begin() + swap.first, P.begin() + swap.second);
				// restart loop
				indices.clear();
			}
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
		return T::Polygon({T::Point(-x, -y), T::Point(-x, y), T::Point(x, y), T::Point(x, -y)});
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
			return T::Polygon({T::Point(-1. * infinity, 0.), T::Point(infinity, 0.), p});
		} else {
			const double scale = 1. / std::sqrt(0.5 * height);
			p.x *= scale;
			p.y *= scale;
			return T::Polygon({T::Point(-0.5 * scale, 0.), T::Point(0.5 * scale, 0.), p});
		}
	}

}
