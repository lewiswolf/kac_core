// core
#include <iostream>	 // cout, endl
#include <stdlib.h>	 // srand
#include <time.h>	 // time

// src
#include "../src/geometry.hpp"

// test
#include "./src/test_modes.cpp"
#include "./src/test_shapes.cpp"

int main() {
	// init random number generator
	testModes();
	testShapes();

	// std::cout << "bessel J_0(4): ";
	// std::cout << bessel(0, 4);
	// std::cout << std::endl << std::endl;

	return 0;
}