// core
#include <iostream>		// cout, endl
#include <stdlib.h>     // srand
#include <time.h>       // time

// src
#include "../src/geometry.cpp"
#include "../src/types.hpp"

int main() {
	// init random number generator
  	srand(time(NULL));

	Vertices v = generateConvexPolygon(7);
	for(int i = 0; i < v.size(); i++) {
       std::cout << "(" << v.at(i).x << ", " << v.at(i).y << ") ";
	}
   	std::cout << std::endl;

	// std::cout << "Hello, World!" << std::endl;

	// Point p1 = Point(1.0, 1.0);
	// std::cout << p1.x << std::endl;
	// std::cout << p1.y << std::endl;

	// Line l = Line(Point(0.1, 0.2), Point(1.2, 1.3));
	// std::cout << l.a.x << ' ' << l.b.y << std::endl;

	// Polygon p = Polygon(v);
   	// std::cout << p.n << std::endl;

	return 0;
}