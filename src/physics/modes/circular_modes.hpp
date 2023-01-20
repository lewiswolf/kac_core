/*
Functions for calculating the linear approximation of the 2-dimensional circular
wave equation.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

// dependencies
#include <boost/math/special_functions/bessel.hpp>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	double besselJ(const int& n, const double& x) {
		/*
		Calculates the bessel function of the first kind J_n(x).
		input:
			n = order of the bessel function
			x = x coordinate
		output:
			y = J_n(x) | y ∈ ℝ
		*/
		return boost::math::cyl_bessel_j(n, x);
	}

	double besselJZero(const double& n, const int& m) {
		/*
		Calculates the mth zero crossing of the bessel functions of the first
		kind.
		input:
			n = order of the bessel function
			m = mth zero
		output:
			z_nm = mth zero crossing of J_n() | z_mn ∈ ℝ
		*/
		return boost::math::cyl_bessel_j_zero(n, m);
	}

	T::Matrix_2D circularAmplitudes(const double& r, const double& theta, const T::Matrix_2D& S) {
		/*
		Calculate the amplitudes of the circular eigenmodes relative to a polar
		strike location.
		input:
			(r, θ) = polar strike location
			S = { z_nm | s ∈ ℝ, J_n(z_nm) = 0, 0 <= n < N, 0 < m <= M }
		output:
			A = {
				abs(J_n(z_nm * r) * (2 ** 0.5) * sin(nθ + π/4))
				| a ∈ ℝ, J_n(z_nm) = 0, 0 <= n < N, 0 < m <= M
			}
		*/

		unsigned int N = S.size();
		unsigned int M = S[0].size();
		T::Matrix_2D A(N, T::Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			double angular = n != 0 ? M_SQRT2 * sin(n * theta + M_PI_4) : 1.0;
			for (unsigned int m = 0; m < M; m++) {
				A[n][m] = abs(besselJ(n, S[n][m] * r) * angular);
			};
		}
		return A;
	}

	T::Matrix_2D circularSeries(const int& N, const int& M) {
		/*
		Calculate the eigenmodes of a circle.
		input:
			N = number of modal orders
			M = number of modes per order
		output:
			S = { z_nm | s ∈ ℝ, J_n(z_nm) = 0, 0 <= n < N, 0 < m <= M }
		*/

		T::Matrix_2D S(N, T::Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) { S[n][m] = besselJZero(n, m + 1); }
		}
		return S;
	}

}