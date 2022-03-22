// core
#include <algorithm>	// shuffle, sort
#include <math.h>       // atan
#include <random>       // default_random_engine
#include <stdlib.h>     // rand, RAND_MAX
#include <time.h>       // time

// src
#include "types.hpp"

Vertices generateConvexPolygon(int const& n) {
	/*
		Generate convex shapes according to Pavel Valtr's 1995 algorithm. Adapted
		from Sander Verdonschot's Java version, found here:
			https://cglab.ca/~sander/misc/ConvexGeneration/ValtrAlgorithm.java
	*/

	// initialise variables
	Vertices vertices;
	double X[n];
	double Y[n];
	double X_rand[n];
	double Y_rand[n];
	unsigned int last_true = 0; 
	unsigned int last_false = 0;
	
	// initialise and sort random coordinates
	for (unsigned int i = 0; i < n; i++) {
		X_rand[i] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		Y_rand[i] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	}
	std::sort(X_rand, X_rand + n);
	std::sort(Y_rand, Y_rand + n);

	// divide the interior points into two chains
	for (unsigned int i = 1; i < n; i++) {
		if (i != n - 1) {
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
	std::shuffle(Y, Y + n, std::default_random_engine(time(NULL)));
	for (unsigned int i = 0; i < n; i++) {
		vertices.push_back(Point(X[i], Y[i]));
	}

	// sort by polar angle
    std::sort(
		vertices.begin(),
		vertices.end(), 
		[&](Point &p1, Point &p2){ return atan2(p1.y, p1.x) < atan2(p2.y, p2.x); }
	);

	// arrange points end to end to form a polygon and center around origin
	double x_min, x_max, y_min, y_max = 0;
	for (unsigned int i = 0; i < n; i++) {
		vertices[i].y += vertices[i].x;
		x_min = fmin(vertices[i].x, x_min);
		x_max = fmax(vertices[i].x, x_max);
		y_min = fmin(vertices[i].y, y_min);
		y_max = fmax(vertices[i].y, y_max);
	}
	for (unsigned int i = 0; i < n; i++) {
		vertices[i].x += ((x_max - x_min) / 2.0) - x_max;
		vertices[i].y += ((y_max - y_min) / 2.0) - y_max;
	}
	return vertices;	  
};