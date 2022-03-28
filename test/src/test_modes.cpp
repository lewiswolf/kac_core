/*
Tests and profiling for /modes.
*/

// core
#include <iostream>
#include <vector>
// src
#include "geometry.hpp"
using namespace geometry;
// test
#include "./utils.hpp"

void testModes() {
	// tests
	std::cout << "Test if besselJ is accurate... ";
	booleanTest(
		std::abs(besselJ(0, 4.2) - -0.37655) < 0.01
		&& std::abs(besselJ(1, 1.2) - 0.498289) < 0.01);
	std::cout << "Test if the 0th mode from calculateLinearModes is f_0... ";
	booleanTest(calculateLinearModes(440.0, 10)[0] == 440.0);
	std::cout << "Test if the 0th element from calculateLinearSeries is 1... ";
	booleanTest(calculateLinearSeries(10)[0] == 1);

	// efficiency
	const int efficiency = 1000;
	std::cout << "\nEfficiency for " << efficiency << " modes...\n";
	{
		Timer timer("	besselJ");
		besselJ(10, 10.2);
	}
	{
		Timer timer("	besselJZero");
		besselJZero(10, 10);
	}
	{
		Timer timer("	calculateCircularModes");
		calculateCircularModes(
			440.0, floor(sqrt(efficiency)), floor(sqrt(efficiency)));
	}
	{
		Timer timer("	calculateCircularSeries");
		calculateCircularSeries(
			floor(sqrt(efficiency)), floor(sqrt(efficiency)));
	}
	{
		Timer timer("	calculateLinearModes");
		calculateLinearModes(440.0, efficiency);
	}
	{
		Timer timer("	calculateLinearSeries");
		calculateLinearSeries(efficiency);
	}
	std::cout << "\n";
}