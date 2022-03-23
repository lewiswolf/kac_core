// core
#include <math.h>	 // atan
#include <stdlib.h>	 // rand, RAND_MAX
#include <time.h>	 // time

#include <algorithm>  // shuffle, sort
#include <random>	  // default_random_engine

Vertices generateConvexPolygon(int const &n) {
	/*
		Generate convex shapes according to Pavel Valtr's 1995 algorithm.
	   Adapted from Sander Verdonschot's Java version, found here:
			https://cglab.ca/~sander/misc/ConvexGeneration/ValtrAlgorithm.java
	*/

	// initialise variables
	Vertices v;
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<double> X_rand;
	std::vector<double> Y_rand;
	unsigned int last_true = 0;
	unsigned int last_false = 0;

	// initialise and sort random coordinates
	for (unsigned int i = 0; i < n; i++) {
		X_rand.push_back(static_cast<double>(rand()) /
						 static_cast<double>(RAND_MAX));
		Y_rand.push_back(static_cast<double>(rand()) /
						 static_cast<double>(RAND_MAX));
		X.push_back(0);
		Y.push_back(0);
	}
	std::sort(X_rand.begin(), X_rand.end());
	std::sort(Y_rand.begin(), Y_rand.end());

	// divide the interior points into two chains
	for (unsigned int i = 1; i < n; i++) {
		if (i != n - 1) {
			// if (rand() % 2 == 1) {
			if (i % 2 == 1) {
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
	std::shuffle(Y.begin(), Y.end(), std::default_random_engine(time(NULL)));
	for (unsigned int i = 0; i < n; i++) {
		v.push_back(Point(X[i], Y[i]));
	}

	// sort by polar angle
	std::sort(v.begin(), v.end(),
			  [](Point p1, Point p2) { return p1.theta < p2.theta; });
	std::cout << std::endl;

	// arrange points end to end to form a polygon
	double x_min, x_max, y_min, y_max = 0;
	double x = 0.0;
	double y = 0.0;
	for (unsigned int i = 0; i < n; i++) {
		Point p = Point(x, y);
		x += v[i].x;
		y += v[i].y;
		v[i] = p;
		x_min = fmin(v[i].x, x_min);
		x_max = fmax(v[i].x, x_max);
		y_min = fmin(v[i].y, y_min);
		y_max = fmax(v[i].y, y_max);
	}
	// center around origin
	for (unsigned int i = 0; i < n; i++) {
		v[i].updateCart(v[i].x + ((x_max - x_min) / 2.0) - x_max,
						v[i].y + ((y_max - y_min) / 2.0) - y_max);
	}
	return v;
};