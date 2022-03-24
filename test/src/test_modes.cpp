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
	std::cout << "Test if the 0th mode from caculateLinearModes is f_0... ";
	booleanTest(modes[0] == 440.0);
	std::cout << "Test if the 0th element from caculateLinearSeries is 1... ";
	booleanTest(calculateLinearSeries(10)[0] == 1);

	// efficiency
	std::cout << "\nEfficiency for 1000 modes...\n";
	std::vector<double> test;
	{
		Timer timer("	calculateLinearModes");
		test = calculateLinearModes(440.0, 1000);
	}
	{
		Timer timer("	calculateLinearSeries");
		test = calculateLinearSeries(1000);
	}
	std::cout << "\n";
}