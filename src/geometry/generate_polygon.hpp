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

static std::default_random_engine random_engine = std::default_random_engine(time(0l));
static std::uniform_int_distribution<long>
	uniform_distribution(std::numeric_limits<long>::min(), std::numeric_limits<long>::max());
static const double rand_max = static_cast<double>(std::numeric_limits<long>::max());

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
			return static_cast<double>(uniform_distribution(random_engine)) / rand_max;
		});
		std::generate(Y_rand.begin(), Y_rand.end(), []() {
			return static_cast<double>(uniform_distribution(random_engine)) / rand_max;
		});
		std::sort(X_rand.begin(), X_rand.end());
		std::sort(Y_rand.begin(), Y_rand.end());
		// divide the interior points into two chains
		for (unsigned long n = 1; n < N; n++) {
			if (n != N - 1) {
				if (rand() % 2 == 1) {
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
		for (unsigned long n = 0; n < N; n++) {
			P[n].x += ((x_max - x_min) / 2.0) - x_max;
			P[n].y += ((y_max - y_min) / 2.0) - y_max;
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
			X.push_back(static_cast<double>(uniform_distribution(random_engine)) / rand_max);
			Y.push_back(static_cast<double>(uniform_distribution(random_engine)) / rand_max);
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
			P.push_back(T::Point(
				static_cast<double>(uniform_distribution(random_engine)) / rand_max,
				static_cast<double>(uniform_distribution(random_engine)) / rand_max
			));
		}
		// 2 opt loop
		std::vector<std::pair<long, long>> indices;
		std::string intersection_type = "";
		bool intersections = true;
		while (intersections) {
		Search_loop:
			for (unsigned long i = 0; i < N - 2; i++) {
				for (unsigned long j = i + 1; j < N; j++) {
					// collect indices of lines which should be crossed
					intersection_type =
						lineIntersection(T::Line(P[i], P[i + 1]), T::Line(P[j], P[(j + 1) % N]))
							.first;
					if (intersection_type == "none" || intersection_type == "vertex") {
						continue;
					} else if (intersection_type == "intersect") {
						indices.push_back(std::make_pair(i + 1, j + 1));
					} else if (intersection_type == "adjacent") {
						indices.push_back(std::make_pair(i, j));
					} else if (intersection_type == "colinear") {
						int min_i = 0 ? P[i].x < P[i + 1].x : 1;
						int max_j = 0 ? P[j].x > P[j + 1].x : 1;
						std::reverse(P.begin() + i + min_i, P.begin() + j + max_j);
						// restart loop
						indices.clear();
						goto Search_loop;
					}
				}
			}
			if (indices.size() > 0) {
				// randomly swap one pair
				std::pair<long, long> swap =
					indices[abs(uniform_distribution(random_engine)) % indices.size()];
				std::reverse(P.begin() + swap.first, P.begin() + swap.second);
				// restart loop
				indices.clear();
				goto Search_loop;
			} else {
				// close loop
				intersections = false;
			}
		}
		return P;
	}

	inline T::Polygon generateUnitRectangle(const double& epsilon) {
		/*
		Define a rectangle with unit area and an aspect ration epsilon.
		*/

		double x = 0.5 * epsilon;
		double y = 0.5 / epsilon;
		return T::Polygon({T::Point(-x, -y), T::Point(-x, y), T::Point(x, y), T::Point(x, -y)});
	}

}