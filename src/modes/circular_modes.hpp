/*
Functions for calculating the linear approximation of the 2-Dimensional circular
wave equation.
*/

#pragma once

// core
#include <vector>
// src
#include "bessel.hpp"

namespace geometry {

	std::vector<std::vector<double>>
	calculateCircularModes(const double& f0, const int& N, const int& M) {
		/*
		Calculate the eigenfrequencies of a circle relative to a given
		fundmental. params: f_0 = fundamental frequency N = number of modes
		*/

		std::vector<std::vector<double>> modes(N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				modes[n][m] = f0 * besselJZero(n, m + 1);
			}
		}
		return modes;
	}

	std::vector<std::vector<double>>
	calculateCircularSeries(const int& N, const int& M) {
		/*
		Calculate the eigenmodes of a circle.
		*/

		std::vector<std::vector<double>> series(N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				series[n][m] = besselJZero(n, m + 1);
			}
		}
		return series;
	}

}