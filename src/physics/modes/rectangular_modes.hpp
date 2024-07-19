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
		const double& x,
		const double& y,
		const unsigned long& N,
		const unsigned long& M,
		const double& epsilon
	) {
		/*
		Calculate the amplitudes of the rectangular eigenmodes relative to a cartesian strike
		location.
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
		for (unsigned long n = 0; n < N; n++) {
			double n_hat = sin((n + 1) * y_hat);
			for (unsigned long m = 0; m < M; m++) { A[n][m] = abs(sin((m + 1) * x_hat) * n_hat); }
		}
		return A;
	}

	inline T::BooleanImage rectangularChladniPattern(
		const double& n,
		const double& m,
		const unsigned long& X,
		const unsigned long& Y,
		const double& tolerance = 0.1
	) {
		/*
		Produce the 2D chladni pattern for a rectangular plate.
		http://paulbourke.net/geometry/chladni/
		input:
			n = nth modal index
			m = mth modal index
			X = length of the X axis
			Y = length of the Y axis
			tolerance = the standard deviation between the calculation and the final pattern
		output:
			M = {
				cos(nπx/X) cos(mπy/Y) - cos(mπx/X) cos(nπy/Y) ≈ 0
			}
		*/

		T::BooleanImage M(X, std::vector<short>(Y, 0));
		for (unsigned long x = 0; x < X; x++) {
			double x_m = cos(m * pi * x / X);
			double x_n = cos(n * pi * x / X);
			for (unsigned long y = 0; y < Y; y++) {
				M[x][y] = abs((x_n * cos(m * pi * y / Y)) - (x_m * cos(n * pi * y / Y))) < tolerance
							? 1
							: 0;
			}
		}
		return M;
	}

	inline T::Matrix_2D
	rectangularSeries(const unsigned long& N, const unsigned long& M, const double& epsilon) {
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
		for (unsigned long n = 0; n < N; n++) {
			double n_hat = pow((n + 1), 2) * epsilon;
			for (unsigned long m = 0; m < M; m++) {
				S[n][m] = sqrt(pow((m + 1), 2) / epsilon + n_hat);
			}
		}
		return S;
	}

}
