/*
Functions for calculating the linear approximation of the 2-dimensional circular
wave equation.
*/

#pragma once

// core
#include <math.h>
#include <numbers>
#include <vector>
using namespace std::numbers;

// dependencies
#include <boost/math/special_functions/bessel.hpp>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	inline double besselJ(const long& n, const double& x) {
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

	inline double besselJZero(const double& n, const long& m) {
		/*
		Calculates the mth zero crossing of the bessel functions of the first kind.
		input:
			n = order of the bessel function
			m = mth zero
		output:
			z_nm = mth zero crossing of J_n() | z_mn ∈ ℝ
		*/

		return boost::math::cyl_bessel_j_zero(n, m);
	}

	inline T::Matrix_2D
	circularAmplitudes(const double& r, const double& theta, const T::Matrix_2D& S) {
		/*
		Calculate the amplitudes of the circular eigenmodes relative to a polar strike location.
		input:
			(r, θ) = polar strike location
			S = { z_nm | s ∈ ℝ, J_n(z_nm) = 0, 0 <= n < N, 0 < m <= M }
		output:
			A = {
				abs(J_n(z_nm * r) * (2 ** 0.5) * sin(nθ + π/4))
				| a ∈ ℝ, J_n(z_nm) = 0, 0 <= n < N, 0 < m <= M
			}
		*/

		const unsigned long N = S.size();
		const unsigned long M = S[0].size();
		const double pi_4 = pi / 4;
		T::Matrix_2D A(N, T::Matrix_1D(M, 0));
		for (unsigned long n = 0; n < N; n++) {
			double angular = n != 0 ? sqrt2 * sin(n * theta + pi_4) : 1.;
			for (unsigned long m = 0; m < M; m++) {
				A[n][m] = abs(boost::math::cyl_bessel_j(n, S[n][m] * r) * angular);
			};
		}
		return A;
	}

	inline T::Matrix_2D circularCymatics(const double& n, const double& m, const unsigned long& H) {
		/*
		Produce the 2D continuos cymatic diagram for a particular mode of a circular domain.
		http://paulbourke.net/geometry/chladni/
		input:
			n = nth modal index
			m = mth modal index
			H = length of the X and Y axis
		output:
			M = {
				J_n(z_nm * r) * (cos(nθ) + sin(nθ)) ≈ 0
			}
		*/

		// interpolate n and z_mn
		double n_round = round(2. * n) / 2.;
		double m_floor = floor(m);
		double z_mn_floor = boost::math::cyl_bessel_j_zero(n, m_floor);
		double z_mn = z_mn_floor
					+ ((boost::math::cyl_bessel_j_zero(n, ceil(m)) - z_mn_floor) * (m - m_floor));
		// calculate pattern
		T::Matrix_2D M(H, T::Matrix_1D(H, 0));
		for (unsigned long x = 0; x < H; x++) {
			double x_prime = (2. * x / H) - 1.;
			for (unsigned long y = 0; y < H; y++) {
				double y_prime = (2. * y / H) - 1.;
				double r = sqrt(pow(x_prime, 2) + pow(y_prime, 2));
				if (r > 1.) {
					continue;
				} else {
					double theta = atan2(y_prime, x_prime);
					M[x][y] = boost::math::cyl_bessel_j(n, z_mn * r)
							* (cos(n_round * theta) + sin(n_round * theta));
				}
			}
		}
		return M;
	}

	inline T::BooleanImage circularChladniPattern(
		const double& n, const double& m, const unsigned long& H, const double& tolerance = 0.1
	) {
		/*
		Produce the 2D Chladni pattern for a circular plate.
		http://paulbourke.net/geometry/chladni/
		input:
			n = nth modal index
			m = mth modal index
			H = length of the X and Y axis
			tolerance = the standard deviation between the calculation and the final pattern
		output:
			B = abs({
				J_n(z_nm * r) * (cos(nθ) + sin(nθ)) ≈ 0
			}) < tolerance
		*/

		// calculate pattern
		T::Matrix_2D M = circularCymatics(n, m, H);
		T::BooleanImage B(H, std::vector<short>(H, 0));
		for (unsigned long x = 0; x < H; x++) {
			for (unsigned long y = 0; y < H; y++) { B[x][y] = abs(M[x][y]) < tolerance ? 1 : 0; }
		}
		return B;
	}

	inline T::Matrix_2D circularSeries(const unsigned long& N, const unsigned long& M) {
		/*
		Calculate the eigenmodes of a circle.
		input:
			N = number of modal orders
			M = number of modes per order
		output:
			S = { z_nm | s ∈ ℝ, J_n(z_nm) = 0, 0 <= n < N, 0 < m <= M }
		*/

		T::Matrix_2D S(N, T::Matrix_1D(M, 0));
		for (unsigned long n = 0; n < N; n++) {
			for (unsigned long m = 0; m < M; m++) {
				S[n][m] = boost::math::cyl_bessel_j_zero((double)n, m + 1);
			}
		}
		return S;
	}

}
