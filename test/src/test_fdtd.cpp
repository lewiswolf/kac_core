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
		p::raisedCosine1D(10, 4, 1.0)[4] == 1.0 && p::raisedCosine2D(10, 10, 4, 4, 1.0)[4][4] == 1.0
	);

	double cfl = 1 / pow(2, 0.5);
	// square simulation
	p::FDTDWaveform2D(
		{{0., 0., 0., 0., 0.},
		 {0., 0., 0., 0., 0.},
		 {0., 0., 0., 0., 0.},
		 {0., 0., 0., 0., 0.},
		 {0., 0., 0., 0., 0.}},
		{{0., 0., 0., 0., 0.},
		 {0., 0., 0., 0., 0.},
		 {0., 0., 1., 0., 0.},
		 {0., 0., 0., 0., 0.},
		 {0., 0., 0., 0., 0.}},
		{{0, 0, 0, 0, 0},
		 {0, 1, 1, 1, 0},
		 {0, 1, 1, 1, 0},
		 {0, 1, 1, 1, 0},
		 {0, 0, 0, 0, 0}},
		pow(cfl, 2),
		2 * (1 - 2 * pow(cfl, 2)),
		1.,
		10,
		{4, 4}
	);
	return 0;
}