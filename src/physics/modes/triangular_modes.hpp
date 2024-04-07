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

	inline T::Matrix_2D equilateralTriangleAmplitudes(
		double u, double v, double w, const unsigned long& N, const unsigned long& M
	) {
		/*
		Calculate the amplitudes of the equilateral triangle eigenmodes relative to a trilinear
		strike location according to Lamé's formula.
		Seth (1940) Transverse Vibrations of Triangular Membranes.
		input:
			( u, v, w ) = trilinear coordinate
			N = number of modal orders
			M = number of modes per order
		output:
			A = {
				abs(sin(nxπ) sin(nyπ) sin(nzπ))
				| a ∈ ℝ, 0 < n <= N, 0 < m <= M
			}
		*/

		u *= pi;
		v *= pi;
		w *= pi;
		T::Matrix_2D A(N, T::Matrix_1D(M, 0));
		for (unsigned long n = 0; n < N; n++) {
			double n_hat = abs(sin((n + 1) * u) * sin((n + 1) * v) * sin((n + 1) * w));
			for (unsigned long m = 0; m < M; m++) { A[n][m] = n_hat; }
		}
		return A;
	}

	inline T::Matrix_2D equilateralTriangleSeries(const unsigned long& N, const unsigned long& M) {
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
		for (unsigned long n = 1; n < N + 1; n++) {
			double n_hat = pow(n, 2);
			for (unsigned long m = 1; m < M + 1; m++) {
				S[n - 1][m - 1] = sqrt(pow(m, 2) + n_hat + (m * n));
			}
		}
		return S;
	}

}
