/*
Functions for calculating the linear approximation of the 2-Dimensional
rectangular wave equation.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

namespace geometry {

	std::vector<std::vector<double>> calculateRectangularAmplitudes(
		const double& x,
		const double& y,
		const int& N,
		const int& M,
		const double& epsilon = 1.0) {
		/*
		Calculate the amplitudes of the eigenmodes relative to rectangular
		strike location.
		*/

		double x_hat = x * M_PI / sqrt(epsilon);
		double y_hat = y * M_PI / sqrt(epsilon);
		std::vector<std::vector<double>> amplitudes(
			N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = sin((n + 1) * y_hat);
			for (unsigned int m = 0; m < M; m++) {
				amplitudes[n][m] = sin((m + 1) * x_hat) * n_hat;
			}
		}
		return amplitudes;
	};

	std::vector<std::vector<double>> calculateRectangularSeries(
		const int& N, const int& M, const double& epsilon = 1.0) {
		/*
		Calculate the eigenmodes of a rectangle.
		*/

		std::vector<std::vector<double>> series(N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = pow((n + 1) * epsilon, 2);
			for (unsigned int m = 0; m < M; m++) {
				series[n][m] = sqrt(pow((m + 1) / epsilon, 2) + n_hat);
			}
		}
		return series;
	};

}