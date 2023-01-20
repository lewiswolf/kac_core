/*
Functions for calculating the linear approximation of the 2-Dimensional
rectangular wave equation.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	T::Matrix_2D equilateralTriangleAmplitudes(
		const double& x, const double& y, const double& z, const int& N, const int& M
	) {
		/*
		Calculate the amplitudes of the equilateral triangle eigenmodes relative to a
		trilinear strike location according to Lamé's formula.
		Seth (1940) Transverse Vibrations of Triangular Membranes.
		input:
			( x, y, z ) = trilinear coordinate
			N = number of modal orders
			M = number of modes per order
		output:
			A = {
				abs(sin(nxπ) sin(nyπ) sin(nzπ))
				| a ∈ ℝ, 0 < n <= N, 0 < m <= M
			}
		*/

		double x_hat = x * M_PI;
		double y_hat = y * M_PI;
		double z_hat = z * M_PI;
		T::Matrix_2D A(N, T::Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = abs(sin((n + 1) * x_hat) * sin((n + 1) * y_hat) * sin((n + 1) * z_hat));
			for (unsigned int m = 0; m < M; m++) { A[n][m] = n_hat; }
		}
		return A;
	};

	T::Matrix_2D equilateralTriangleSeries(const int& N, const int& M) {
		/*
		Calculate the eigenmodes of an equilateral triangle according to Lamé's formula.
		Seth (1940) Transverse Vibrations of Triangular Membranes.
		input:
			N = number of modal orders
			M = number of modes per order
		output:
			S = {
				(m ** 2 + n ** 2 + mn) ** 0.5
				| s ∈ ℝ, 0 < n <= N, 0 < m <= M
			}
		*/

		T::Matrix_2D S(N, T::Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double n_hat = pow((n + 1), 2);
			for (unsigned int m = 0; m < M; m++) {
				S[n][m] = sqrt(pow((m + 1), 2) + n_hat + m * n);
			}
		}
		return S;
	};

}