#pragma once

// core
#include <algorithm>
#include <vector>

std::vector<double> calculateLinearModes(const double& f0, const int& N) {
	/*
	Calculate the harmonic series of a given fundamental.

	params:
		f0	-	fundamental frequency
		N	-	number of modes
	*/

	std::vector<double> modes(N);
	for (unsigned int i = 0; i < N; i++) { modes[i] = f0 * (i + 1); };
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