#pragma once

// core
#include <cmath>
#include <vector>

namespace geometry {

	std::vector<std::vector<double>> calculateRectangularModes(
		const int& N,
		const int& M,
		const double& f_0,
		const double& epsilon = 1.0) {
		/*
		Calculate the eigenmodes of a rectangle.
		*/

		std::vector<std::vector<double>> modes(N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				modes[n][m] =
					f_0 * sqrt(pow((m + 1) / epsilon, 2) + pow((n + 1), 2));
			}
		}
		return modes;
	}

	std::vector<std::vector<double>> calculateRectangularSeries(
		const int& N, const int& M, const double& epsilon = 1.0) {
		/*
		Calculate the eigenmodes of a rectangle.
		*/

		std::vector<std::vector<double>> series(N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				series[n][m] =
					sqrt(pow((m + 1) / epsilon, 2) + pow((n + 1), 2));
			}
		}
		return series;
	}

}