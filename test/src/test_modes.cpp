/*
Tests and profiling for /modes.
*/

// core
#include <iostream>	   // cout
#include <stdlib.h>	   // rand, RAND_MAX
#include <vector>
// src
#include "geometry.hpp"
using namespace geometry;
// test
#include "./utils.hpp"

void testModes() {
	// tests
	std::cout << "besselJ is accurate... ";
	booleanTest(
		std::abs(besselJ(0, 4.2) - -0.37655) < 0.01
		&& std::abs(besselJ(1, 1.2) - 0.498289) < 0.01
	);
	std::cout << "the 0th mode from calculateLinearModes is f_0... ";
	booleanTest(calculateLinearModes(440.0, 10)[0] == 440.0);
	std::cout << "the 0th element from calculateLinearSeries is 1... ";
	booleanTest(calculateLinearSeries(10)[0] == 1);

	// efficiency
	const int efficiency = 1000;
	std::cout << "\nEfficiency for " << efficiency << " modes...\n";
	{
		Timer timer("	besselJ");
		for (unsigned int i = 0; i < floor(sqrt(efficiency)); i++) {
			for (unsigned int j = 0; j < floor(sqrt(efficiency)); j++) {
				besselJ(i, j);
			}
		}
	}
	{
		Timer timer("	besselJZero");
		for (unsigned int i = 0; i < floor(sqrt(efficiency)); i++) {
			for (unsigned int j = 0; j < floor(sqrt(efficiency)); j++) {
				besselJZero(i, j);
			}
		}
	}
	{
		Timer timer("	calculateCircularModes");
		calculateCircularModes(
			440.0, floor(sqrt(efficiency)), floor(sqrt(efficiency))
		);
	}
	{
		Timer timer("	calculateCircularSeries");
		calculateCircularSeries(
			floor(sqrt(efficiency)), floor(sqrt(efficiency))
		);
	}
	{
		Timer timer("	calculateLinearAmplitudes");
		calculateLinearAmplitudes(
			static_cast<double>(rand()) / static_cast<double>(RAND_MAX),
			efficiency
		);
	}
	{
		Timer timer("	calculateLinearModes");
		calculateLinearModes(440.0, efficiency);
	}
	{
		Timer timer("	calculateLinearSeries");
		calculateLinearSeries(efficiency);
	}
	{
		Timer timer("	calculateRectangularAmplitudes");
		calculateRectangularAmplitudes(
			0.5, 0.5, floor(sqrt(efficiency)), floor(sqrt(efficiency))
		);
	}
	{
		Timer timer("	calculateRectangularSeries");
		calculateRectangularSeries(
			floor(sqrt(efficiency)), floor(sqrt(efficiency))
		);
	}
	std::cout << "\n";
}