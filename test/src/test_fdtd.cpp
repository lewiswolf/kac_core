/*
Tests for /fdtd.
*/

// src
#include <kac_core.hpp>
namespace p = kac_core::physics;

// test
#include "./utils.hpp"

int main() {
	booleanTest(
		"raisedCosine is accurate",
		p::raisedCosine1D(10, 4, 1.0)[4] == 1.0
			&& p::raisedCosine2D(10, 10, 4, 4, 1.0)[4][4] == 1.0
	);
	return 0;
}