// core
#pragma once
#include <array>
#include <vector>

typedef struct Point {
	/*
		A cartesian product.
	*/
	double x = 0.0;
	double y = 0.0;

	Point() {};
	Point(double x, double y): x(x), y(y) {};
} Point;

typedef struct Line {
	/*
		A line from point a to point b.
	*/
	Point a;
	Point b;

	Line() {};
	Line(Point a, Point b): a(a), b(b) {};
} Line;

// An array of Points
struct Vertices : public std::vector<Point> {
	/*
		An array of points.
	*/
	std::vector<std::array<double, 2>> convertVerticesToVector() {
		// covert to vector of arrays
		std::vector<std::array<double, 2>> out;
		for (int i = 0; i < size(); i++) {
			out.push_back({{at(i).x, at(i).y}});
		}
		return out;
	}
};

// typedef struct Circle {

// } Circle;

// typedef struct Polygon {
// 	/*
// 		An ordered set of n vertices.
// 	*/
// 	int n;				// number of vertices
// 	Vertices vertices;	// dynamic array of vertices

// 	Polygon() {};
// 	Polygon(Vertices v): vertices(v) {
// 		n = v.size();
// 	};
// } Polygon;