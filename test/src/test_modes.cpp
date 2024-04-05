/*
Tests for /modes.
*/

// core
#include <math.h>

// src
#include <kac_core.hpp>
namespace p = kac_core::physics;

// test
#include "./utils.hpp"

int main() {
	/*
	Test Bessel function.
	*/
	booleanTest(
		"besselJ is accurate",
		abs(p::besselJ(0, 4.2) - -0.37655) < 0.01 && abs(p::besselJ(1, 1.2) - 0.498289) < 0.01
	);

	/*
	Test linear modes.
	*/
	booleanTest("the 0th element from linearSeries is 1", p::linearSeries(10)[0] == 1);

	return 0;
}
