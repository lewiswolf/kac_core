// core
#include <iostream>
#include <vector>

// src
#include "../../src/geometry.hpp"

// test
#include "./utils.hpp"

void testModes() {
	// tests
	std::vector<double> modes = calculateLinearModes(440.0, 10);
	std::cout << "Test if besselJ is accurate... ";
	booleanTest(std::abs(besselJ(0, 4.2) - -0.37655) < 0.01 &&
				std::abs(besselJ(1, 1.2) - 0.498289) < 0.01);
	std::cout << "Test if the 0th mode from caculateLinearModes is f_0... ";
	booleanTest(modes[0] == 440.0);
	std::cout << "Test if the 0th element from caculateLinearSeries is 1... ";
	booleanTest(calculateLinearSeries(10)[0] == 1);

	// efficiency
	const int efficiency = 1000;
	std::cout << "\nEfficiency for " << efficiency << " modes...\n";
	std::vector<double> test;
	{
		Timer timer("	besselJ");
		besselJ(10, 10.2);
	}
	{
		Timer timer("	calculateLinearModes");
		test = calculateLinearModes(440.0, efficiency);
	}
	{
		Timer timer("	calculateLinearSeries");
		test = calculateLinearSeries(efficiency);
	}
	std::cout << "\n";
}