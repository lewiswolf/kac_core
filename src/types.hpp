#pragma once

// core
#include <math.h>
#include <vector>

typedef struct Point {
	/*
	A point on the Euclidian plane.
	*/

	// vars
	double x = 0.0;
	double y = 0.0;
	double r = 0.0;
	double theta = 0.0;

	// constructors
	Point(){};
	Point(double x, double y) { updateCart(x, y); };

	// methods
	void updateCart(double new_x, double new_y) {
		/*
		Update the point using cartesian coordinates.
		*/

		x = new_x;
		y = new_y;
		r = pow(pow(x, 2) + pow(y, 2), 0.5);
		theta = atan2(y, x);
	}

	void updatePol(double new_r, double new_theta) {
		/*
		Update the point using polar coordinates.
		*/

		r = new_r;
		theta = new_theta;
		x = r * cos(theta);
		y = r * sin(theta);
	}
} Point;

typedef struct Line {
	/*
	A straight line from point a to point b.
	*/

	// vars
	Point a;
	Point b;

	// constructors
	Line(){};
	Line(Point a, Point b): a(a), b(b){};
} Line;

typedef std::vector<Point> Vertices;

typedef struct Polygon {
} Polygon;

typedef struct Circle {
	/*
	A circle defined on the euclidian plane.
	*/

	// vars
	double r = 1.0;	 // radius
	Point origin;	 // center

	// constructors
	Circle(){};
	Circle(double r): r(r){};
	Circle(Point origin): origin(origin){};
	Circle(double r, Point origin): r(r), origin(origin){};
} Circle;

typedef struct Ellipse {
} Ellipse;