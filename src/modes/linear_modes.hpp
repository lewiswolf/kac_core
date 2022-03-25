/*
Functions for calculating the linear approximation of a 1-Dimensional wave
equation.
*/

#pragma once
// core
#include <vector>

std::vector<double> calculateLinearModes(const double& f_0, const int& N) {
	/*
	Calculate the harmonic series of a given fundmental.
	params:
		f_0 = fundamental frequency
		N = number of modes
	*/

	std::vector<double> modes(N);
	for (unsigned int i = 0; i < N; i++) { modes[i] = f_0 * (i + 1); };
	return modes;
}

std::vector<double> calculateLinearSeries(const int& N) {
	/*
	Calculate the relationship between the modes of the harmonic series.
	params:
		N	-	number of modes
	*/

	std::vector<double> series(N);
	for (unsigned int i = 0; i < N; i++) { series[i] = i + 1; };
	return series;
}