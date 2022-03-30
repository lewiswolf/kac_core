/*
Functions for calculating the linear approximation of a 1-Dimensional wave
equation.
*/

#pragma once

// core
#include <vector>

namespace geometry {

	std::vector<double> calculateLinearModes(const double& f_0, const int& N) {
		/*
		Calculate the harmonic series of a given fundmental.
		params:
			f_0 = fundamental frequency
			N = number of modes
		output:
			F = { (f_0 * n) | n < N }
		*/

		std::vector<double> F(N);
		for (unsigned int n = 0; n < N; n++) { F[n] = f_0 * (n + 1); };
		return F;
	}

	std::vector<double> calculateLinearSeries(const int& N) {
		/*
		Calculate the relationship between the modes of the harmonic series.
		params:
			N = number of modes
		output:
			S = { n | n < N }
		*/

		std::vector<double> S(N);
		for (unsigned int n = 0; n < N; n++) { S[n] = n + 1; };
		return S;
	}

}