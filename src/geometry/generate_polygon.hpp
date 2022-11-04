/*
Functions for generating polygons.
*/

#pragma once

// core
#include <algorithm>	// generate, min, max, shuffle, sort
#include <random>		// default_random_engine
#include <stdlib.h>		// rand, RAND_MAX
#include <time.h>		// time

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	T::Vertices generateConvexPolygon(const int& N, time_t seed = time(NULL)) {
		/*
		Generate convex shapes according to Pavel Valtr's 1995 algorithm.
		Adapted from Sander Verdonschot's Java version, found here:
		https://cglab.ca/~sander/misc/ConvexGeneration/ValtrAlgorithm.java
		input:
			N = the number of vertices
			seed? = the seed for the random number generators
		output:
			V = a polygon of N random vertices
		*/

		// initialise variables
		T::Vertices V;
		std::vector<double> X(N);
		std::vector<double> Y(N);
		std::vector<double> X_rand(N);
		std::vector<double> Y_rand(N);
		unsigned int last_true = 0;
		unsigned int last_false = 0;
		// initialise and sort random coordinates
		srand(seed);
		std::generate(X_rand.begin(), X_rand.end(), []() {
			return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		});
		std::generate(Y_rand.begin(), Y_rand.end(), []() {
			return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		});
		std::sort(X_rand.begin(), X_rand.end());
		std::sort(Y_rand.begin(), Y_rand.end());
		// divide the interior points into two chains
		for (unsigned int i = 1; i < N; i++) {
			if (i != N - 1) {
				if (rand() % 2 == 1) {
					X[i] = X_rand[i] - X_rand[last_true];
					Y[i] = Y_rand[i] - Y_rand[last_true];
					last_true = i;
				} else {
					X[i] = X_rand[last_false] - X_rand[i];
					Y[i] = Y_rand[last_false] - Y_rand[i];
					last_false = i;
				}
			} else {
				X[0] = X_rand[i] - X_rand[last_true];
				Y[0] = Y_rand[i] - Y_rand[last_true];
				X[i] = X_rand[last_false] - X_rand[i];
				Y[i] = Y_rand[last_false] - Y_rand[i];
			}
		}
		// randomly combine x and y
		shuffle(Y.begin(), Y.end(), std::default_random_engine(seed));
		for (unsigned int i = 0; i < N; i++) { V.push_back(T::Point(X[i], Y[i])); }
		// sort by polar angle
		sort(V.begin(), V.end(), [](T::Point& p1, T::Point& p2) {
			return p1.theta() < p2.theta();
		});
		// arrange points end to end to form a polygon
		double x_min, x_max, y_min, y_max = 0;
		double x = 0.0;
		double y = 0.0;
		for (unsigned int i = 0; i < N; i++) {
			T::Point p = T::Point(x, y);
			x += V[i].x;
			y += V[i].y;
			V[i] = p;
			x_min = std::min(V[i].x, x_min);
			x_max = std::max(V[i].x, x_max);
			y_min = std::min(V[i].y, y_min);
			y_max = std::max(V[i].y, y_max);
		}
		// center around origin
		for (unsigned int i = 0; i < N; i++) {
			V[i].x += ((x_max - x_min) / 2.0) - x_max;
			V[i].y += ((y_max - y_min) / 2.0) - y_max;
		}
		return V;
	}

}