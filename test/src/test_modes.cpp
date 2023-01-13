/*
Tests for /modes.
*/

// src
#include <kac_core.hpp>
namespace p = kac_core::physics;

// test
#include "./utils.hpp"

int main() {
	booleanTest(
		"besselJ is accurate",
		std::abs(p::besselJ(0, 4.2) - -0.37655) < 0.01
			&& std::abs(p::besselJ(1, 1.2) - 0.498289) < 0.01
	);
	booleanTest(
		"the 0th element from calculateLinearSeries is 1", p::calculateLinearSeries(10)[0] == 1
	);
	return 0;
}