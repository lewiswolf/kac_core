/*
Functions for calculating the linear approximation of the 2-Dimensional
rectangular wave equation.
*/

#pragma once

// core
#include <math.h>
#include <numbers>
#include <vector>
using namespace std::numbers;

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	inline T::Matrix_2D rectangularAmplitudes(
		const double& x, const double& y, const int& N, const int& M, const double& epsilon
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
				abs(sin(mxπ / (Є ** 0.5)) sin(nyπ * (Є ** 0.5)))
				| a ∈ ℝ, 0 < n <= N, 0 < m <= M
			}
		*/

		double x_hat = x * pi / sqrt(epsilon);
		double y_hat = y * pi * sqrt(epsilon);
		T::Matrix_2D A(N, T::Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = sin((n + 1) * y_hat);
			for (unsigned int m = 0; m < M; m++) { A[n][m] = abs(sin((m + 1) * x_hat) * n_hat); }
		}
		return A;
	}

	inline T::Matrix_2D rectangularSeries(const int& N, const int& M, const double& epsilon) {
		/*
		Calculate the eigenmodes of a rectangle.
		input:
			N = number of modal orders
			M = number of modes per order
			epsilon = aspect ratio of the rectangle
		output:
			S = {
				((m ** 2 / Є) + (Єn ** 2)) ** 0.5
				| s ∈ ℝ, 0 < n <= N, 0 < m <= M
			}
		*/

		T::Matrix_2D S(N, T::Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = pow((n + 1), 2) * epsilon;
			for (unsigned int m = 0; m < M; m++) {
				S[n][m] = sqrt(pow((m + 1), 2) / epsilon + n_hat);
			}
		}
		return S;
	}

}
