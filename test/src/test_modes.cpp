/*
Tests for /modes.
*/

// core
#include <cmath>

// src
#include <kac_core.hpp>
namespace p = kac_core::physics;

// test
#include "./utils.hpp"

int main() {
	/*
	Test linear modes.
	*/
	booleanTest("the 0th element from linearSeries is 1", p::linearSeries(10)[0] == 1);

	return 0;
}
