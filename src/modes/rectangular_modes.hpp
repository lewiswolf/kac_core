/*
Functions for calculating the linear approximation of the 2-Dimensional
rectangular wave equation.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
// src
#include "../types.hpp"

namespace geometry {

	Matrix_2D calculateRectangularAmplitudes(
		const double& x,
		const double& y,
		const int& N,
		const int& M,
		const double& epsilon = 1.0
	) {
		/*
		Calculate the amplitudes of the rectangular eigenmodes relative to a
		cartesian strike location.
		input:
			( x , y ) = cartesian product
			N = number of modal orders
			M = number of modes per order
			epsilon = aspect ratio of the rectangle
		output:
			A = {
				(sin(nyπ / (epsilon ** 0.5) sin(mxπ / (epsilon ** 0.5))
				| a ∈ ℝ, 0 < n <= N, 0 < m <= M
			}
		*/

		double x_hat = x * M_PI / sqrt(epsilon);
		double y_hat = y * M_PI / sqrt(epsilon);
		Matrix_2D A(N, Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = sin((n + 1) * y_hat);
			for (unsigned int m = 0; m < M; m++) {
				A[n][m] = sin((m + 1) * x_hat) * n_hat;
			}
		}
		return A;
	};

	Matrix_2D calculateRectangularSeries(
		const int& N, const int& M, const double& epsilon = 1.0
	) {
		/*
		Calculate the eigenmodes of a rectangle.
		input:
			N = number of modal orders
			M = number of modes per order
			epsilon = aspect ratio of the rectangle
		output:
			S = {
				((m / epsilon)^2 + (n * epsilon)^2) ** 0.5
				| s ∈ ℝ, 0 < n <= N, 0 < m <= M
			}
		*/

		Matrix_2D S(N, Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = pow((n + 1) * epsilon, 2);
			for (unsigned int m = 0; m < M; m++) {
				S[n][m] = sqrt(pow((m + 1) / epsilon, 2) + n_hat);
			}
		}
		return S;
	};

}