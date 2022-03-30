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
	calculateCircularModes(const double& f_0, const int& N, const int& M) {
		/*
		Calculate the eigenfrequencies of a circle relative to a given
		fundmental.
		params:
			f_0 = fundamental frequency
			N = number of modal order
			M = number of modes per order
		output:
			F = { (f_0 * z_nm) ∈ ℝ | J_n(z_nm) = 0, n < N, m < M }
		*/

		std::vector<std::vector<double>> F(N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				F[n][m] = f_0 * besselJZero(n, m + 1);
			}
		}
		return F;
	}

	std::vector<std::vector<double>>
	calculateCircularSeries(const int& N, const int& M) {
		/*
		Calculate the eigenmodes of a circle.
		params:
			N = number of modal orders
			M = number of modes per order
		output:
			S = { z_nm ∈ ℝ | J_n(z_nm) = 0, n < N, m < M }
		*/

		std::vector<std::vector<double>> S(N, std::vector<double>(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				S[n][m] = besselJZero(n, m + 1);
			}
		}
		return S;
	}

}