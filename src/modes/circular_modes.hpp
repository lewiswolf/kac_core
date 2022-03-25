/*
Functions for calculating the linear approximation of the 2-Dimensional circular
wave equation.
*/

#pragma once
// core
#include <vector>

std::vector<std::vector<double>>
calculateCircularModes(double const& f0, int const& N, int const& M) {
	std::vector<std::vector<double>> modes(N, std::vector<double>(M, 0));
	/*
	Calculate the eigenfrequencies of a circle relative to a given fundmental.
	params:
		f_0 = fundamental frequency
		N = number of modes
	*/

	for (unsigned int n = 0; n < N; n++) {
		for (unsigned int m = 0; m < M; m++) {
			modes[n][m] = f0 * besselJZero(n, m + 1);
		}
	}
	return modes;
}

std::vector<std::vector<double>>
calculateCircularSeries(int const& N, int const& M) {
	std::vector<std::vector<double>> series(N, std::vector<double>(M, 0));
	/*
	Calculate the eigenmodes of a circle.
	*/

	for (unsigned int n = 0; n < N; n++) {
		for (unsigned int m = 0; m < M; m++) {
			series[n][m] = besselJZero(n, m + 1);
		}
	}
	return series;
}