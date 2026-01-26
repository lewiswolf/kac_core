/*
Functions for calculating the linear approximation of the 2-dimensional circular
wave equation.
*/

#pragma once

// core
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

// dependencies
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/tools/roots.hpp>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

// root finding constants
const auto root_finder_tolerance = boost::math::tools::eps_tolerance<double>(50);

namespace kac_core::physics {

	inline T::Matrix_2D
	circularAmplitudes(const double& r, const double& theta, const T::Matrix_2D& S) {
		/*
		Calculate the spatial eigenfunction of a circular 2-dimensional domain relative to a
		polar excitation. The boundary conditions for this spatial eigenfunction are
		determined by the input series of wavenumbers λ_mn.
		input:
			(r, θ) = polar excitation
			S = { λ_mn | λ ∈ ℝ }
		output:
			α_mn = {
				J_m(λ_mn * r√π) * √2 * sin(mθ + π/4)
				| α ∈ ℝ, m ∈ [0, M), n ∈ (0, N]
			}
		*/

		const std::size_t M = S.size();
		const std::size_t N = S[0].size();
		T::Matrix_2D A(M, T::Matrix_1D(N, 0.));
		const double pi_4 = std::numbers::pi * 0.25;
		for (std::size_t m = 0; m < M; m++) {
			double angular = m == 0 ? 1. : std::numbers::sqrt2 * std::sin(m * theta + pi_4);
			for (std::size_t n = 0; n < N; n++) {
				A[m][n] = boost::math::cyl_bessel_j(m, S[m][n] * r * std::sqrt(std::numbers::pi))
						* angular;
			};
		}
		return A;
	}

	inline T::Matrix_2D circularCymatics(
		const double& m,
		const double& n,
		const std::size_t& H,
		const bool boundary_conditions = true
	) {
		/*
		Produce the cymatic diagram of a 2-dimensional circular domain for a particular mode λ_mn.
		For creative use, m and n have been defined as real numbers to create continuous animations,
		however for analytics these should be interpreted as integers.
		http://paulbourke.net/geometry/chladni/
		input:
			m = mth modal index
			n = nth modal index
			H = length of the X and Y axes
			boundary_conditions = (true = fixed, false = free)
		output:
			U_rθ = {
				J_n(z_nm * r) * (cos(nθ) + sin(nθ))
				| U ∈ ℝ^2
			}
		*/

		// interpolate n and z_mn
		double z_mn = 0.;
		if (boundary_conditions) {
			const double n_hat = n + 1.;
			const double n_floor = std::floor(n_hat);
			const double z_mn_floor = boost::math::cyl_bessel_j_zero(m, n_floor);
			z_mn = z_mn_floor
				 + ((boost::math::cyl_bessel_j_zero(m, std::ceil(n_hat)) - z_mn_floor)
					* (n_hat - n_floor));
		} else {
			const double n_floor = std::floor(n);
			const auto j_prime = [m](double _x) { return boost::math::cyl_bessel_j_prime(m, _x); };
			const double lower_bound = n_floor == 0. ? std::numeric_limits<double>::epsilon()
													 : boost::math::cyl_bessel_j_zero(m, n_floor);
			const double mid_bound = boost::math::cyl_bessel_j_zero(m, n_floor + 1);
			const double upper_bound = boost::math::cyl_bessel_j_zero(m, n_floor + 2);
			const double z_mn_floor = (m == 0. && n < 1.)
										? 0.
										: boost::math::tools::bisect(
											  j_prime, lower_bound, mid_bound, root_finder_tolerance
										  )
											  .first;
			const double z_mn_ceil =
				boost::math::tools::bisect(j_prime, mid_bound, upper_bound, root_finder_tolerance)
					.first;
			z_mn = z_mn_floor + ((z_mn_ceil - z_mn_floor) * (n - n_floor));
		}
		// calculate pattern
		T::Matrix_2D U(H, T::Matrix_1D(H, 0.));
		const double H_2 = 2. / H;
		const double m_round = std::round(2. * m) * 0.5;
		for (std::size_t x = 0; x < H; x++) {
			double x_prime = (x * H_2) - 1.;
			for (std::size_t y = 0; y < H; y++) {
				double y_prime = (y * H_2) - 1.;
				double r = std::hypot(x_prime, y_prime);
				if (r <= 1.) {
					double theta = std::atan2(y_prime, x_prime);
					U[x][y] = boost::math::cyl_bessel_j(m, z_mn * r)
							* (std::cos(m_round * theta) + std::sin(m_round * theta));
				}
			}
		}
		return U;
	}

	inline T::Matrix_2D circularSeries(
		const std::size_t& M, const std::size_t& N, const bool boundary_conditions = true
	) {
		/*
		Calculate the wavenumbers of a 2-dimensional circular domain.
		input:
			M = number of modes across the Mth axis
			N = number of modes across the Nth axis
			boundary_conditions = (true = fixed, false = free)
		output:
			z_mn = {
				J_m(z_mn) = 0 					dirichlet boundary condition
				J'_m(z_mn) = 0 					neumann boundary condition
				| m ∈ [0, M), n ∈ (0, N]
			}
			λ_mn { z_mn / √π | λ ∈ ℝ }
		*/

		T::Matrix_2D S(M, T::Matrix_1D(N, 0.));
		if (boundary_conditions) {
			for (std::size_t m = 0; m < M; m++) {
				double nu = static_cast<double>(m);
				for (std::size_t n = 0; n < N; n++) {
					S[m][n] =
						boost::math::cyl_bessel_j_zero(nu, n + 1) / std::sqrt(std::numbers::pi);
				}
			}
		} else {
			double lower_bound = 0., upper_bound = 0.;
			for (std::size_t m = 0; m < M; m++) {
				double nu = static_cast<double>(m);
				auto j_prime = [nu](double _x) { return boost::math::cyl_bessel_j_prime(nu, _x); };
				for (std::size_t n = 0; n < N; n++) {
					lower_bound = n == 0 ? std::numeric_limits<double>::epsilon() : upper_bound;
					upper_bound = boost::math::cyl_bessel_j_zero(nu, n + 1);
					// rigid body mode
					if (m == 0 && n == 0) {
						continue;
					}
					S[m][n] = boost::math::tools::bisect(
								  j_prime, lower_bound, upper_bound, root_finder_tolerance
							  )
								  .first
							/ std::sqrt(std::numbers::pi);
				}
			}
		}
		return S;
	}

}
