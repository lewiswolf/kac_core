/*
Tests for /fdtd.
*/

// src
#include <kac_core.hpp>
namespace p = kac_core::physics;

// test
#include "./utils.hpp"

int main() {
	/*
	Test raisedCosine.
	*/
	booleanTest(
		"raisedCosine is accurate",
		p::raisedCosine1D(0.5, 0.1, 11)[5] == 1.
			&& p::raisedCosine2D(T::Point(0.5, 0.5), 0.1, 11, 11)[5][5] == 1.
	);

	/*
	Test raisedTriangle.
	*/
	booleanTest(
		"raisedTriangle is accurate",
		p::raisedTriangle1D(0.5, 0.1, 0.1, 11)[5] == 1.
			&& p::raisedTriangle2D(T::Point(0.5, 0.5), 0.1, 0.1, 0.1, 0.1, 11, 11)[5][5] == 1.
	);

	/*
	Create a square FDTD simulation.
	*/
	double cfl_2 = pow(1 / pow(2, 0.5), 2.);
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
		{{0, 0, 0, 0, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 1, 1, 1, 0}, {0, 0, 0, 0, 0}},
		cfl_2,
		2 - 4 * cfl_2,
		1.,
		10,
		T::Point(0.5, 0.5)
	);

	return 0;
}
