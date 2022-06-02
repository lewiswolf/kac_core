/*
Tests for /modes.
*/

// src
#include <kac_core.hpp>
using namespace geometry;

// test
#include "./utils.hpp"

int main() {
	booleanTest(
		"besselJ is accurate",
		std::abs(besselJ(0, 4.2) - -0.37655) < 0.01
			&& std::abs(besselJ(1, 1.2) - 0.498289) < 0.01
	);
	booleanTest(
		"the 0th mode from calculateLinearModes is f_0",
		calculateLinearModes(440.0, 10)[0] == 440.0
	);
	booleanTest(
		"the 0th element from calculateLinearSeries is 1",
		calculateLinearSeries(10)[0] == 1
	);
	return 0;
}