/*
Tests for /fdtd.
*/

// core
#include <algorithm>

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
	double cfl_1 = 1.;
	T::Matrix_1D waveform = p::FDTDWaveform1D(
		{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
		{0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.},
		cfl_1,
		2 - 2 * cfl_1,
		1.,
		1000,
		0.5
	);
	auto [min_waveform, max_waveform] = std::minmax_element(begin(waveform), end(waveform));
	booleanTest("FDTDWaveform1D does not explode", *min_waveform >= -1. && *max_waveform <= 1.);

	/*
	Create a square FDTD simulation.
	*/
	double cfl_2 = pow(1 / pow(2, 0.5), 2.);
	waveform = p::FDTDWaveform2D(
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
		1000,
		T::Point(0.5, 0.5)
	);
	auto [min_waveform2, max_waveform2] = std::minmax_element(begin(waveform), end(waveform));
	booleanTest("FDTDWaveform2D does not explode", *min_waveform2 >= -1. && *max_waveform2 <= 1.);

	return 0;
}
