/*
Functions for calculating the linear approximation of the 2-Dimensional circular
wave equation.
*/

#pragma once

// core
#include <vector>

// dependencies
#include <boost/math/special_functions/bessel.hpp>

// src
#include "../../types.hpp"
using namespace kac_core::types;

namespace kac_core { namespace physics {

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

	Matrix_2D
	calculateCircularModes(const double& f_0, const int& N, const int& M) {
		/*
		Calculate the eigenfrequencies of a circle relative to a given
		fundamental.
		input:
			f_0 = fundamental frequency
			N = number of modal order
			M = number of modes per order
		output:
			F = { (f_0 * z_nm) | f ∈ ℝ, J_n(z_nm) = 0, n < N, 0 < m <= M }
		*/

		Matrix_2D F(N, Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				F[n][m] = f_0 * besselJZero(n, m + 1);
			}
		}
		return F;
	}

	Matrix_2D calculateCircularSeries(const int& N, const int& M) {
		/*
		Calculate the eigenmodes of a circle.
		input:
			N = number of modal orders
			M = number of modes per order
		output:
			S = { z_nm | s ∈ ℝ, J_n(z_nm) = 0, n < N, 0 < m <= M }
		*/

		Matrix_2D S(N, Matrix_1D(M, 0));
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				S[n][m] = besselJZero(n, m + 1);
			}
		}
		return S;
	}

}}