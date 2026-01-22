/*
Functions for calculating the linear approximation of the 2-dimensional rectangular wave equation.
*/

#pragma once

// core
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	inline T::Matrix_2D rectangularAmplitudes(
		const double& x,
		const double& y,
		const std::size_t& M,
		const std::size_t& N,
		const double& epsilon,
		const std::array<bool, 4> boundary_conditions = {true, true, true, true}
	) {
		/*
		Calculate the spatial eigenfunction of a rectangular 2-dimensional domain relative to a
		cartesian strike location.
		input:
			(x, y) = cartesian strike location
			M = number of modes across the Mth axis
			N = number of modes across the Nth axis
			epsilon = aspect ratio of the rectangle
			boundary_conditions = boolean array indicating the boundary conditions
				(true = fixed, false = free)
				1: x-axis minima boundary condition
				2: x-axis maxima boundary condition
				3: y-axis minima boundary condition
				4: y-axis maxima boundary condition
		output:
			X_m = {
				sin((m + 1)xπ / √Є),	dirichlet boundary condition
				cos(mxπ / √Є),			neumann boundary condition
				sin((m + 0.5)xπ / √Є),	mixed boundary conditions
				| m ∈ [0, M)
			}
			Y_n = {
				sin((n + 1)yπ√Є),		dirichlet boundary condition
				cos(nyπ√Є),				neumann boundary condition
				sin((n + 0.5)yπ√Є),		mixed boundary conditions
				| n ∈ [0, N)
			}
			α_mn = { X_m * Y_n | α ∈ ℝ }
		*/

		const double epsilon_root = std::sqrt(epsilon);
		const double x_hat = x * std::numbers::pi / epsilon_root;
		const double y_hat = y * std::numbers::pi * epsilon_root;
		T::Matrix_2D A(M, T::Matrix_1D(N, 0.));
		// boundary condition lambda
		auto BCLambda = [](const std::size_t i,
						   const double scalar,
						   const bool minima,
						   const bool maxima) -> double {
			if (minima && maxima) {
				// dirichlet boundary
				return std::sin((i + 1.) * scalar);
			} else if (!minima && !maxima) {
				// neumann boundary
				return std::cos(i * scalar);
			} else {
				// mixed boundary
				return std::sin((i + 0.5) * scalar);
			}
		};
		// produce the spatial eigenfunction
		for (std::size_t m = 0; m < M; m++) {
			double X_m = BCLambda(m, x_hat, boundary_conditions[0], boundary_conditions[1]);
			for (std::size_t n = 0; n < N; n++) {
				A[m][n] = X_m * BCLambda(n, y_hat, boundary_conditions[2], boundary_conditions[3]);
			}
		}
		return A;
	}

	inline T::Matrix_2D rectangularCymatics(
		const double& m,
		const double& n,
		const std::size_t& X,
		const std::size_t& Y,
		const std::array<bool, 4> boundary_conditions = {true, true, true, true}
	) {
		/*
		Produce the cymatic diagram of a 2-dimensional rectangular domain for a particular
		mode λ_mn.
		input:
			m = mth modal index
			n = nth modal index
			X = length of the X axis
			Y = length of the Y axis
			boundary_conditions = boolean array indicating the boundary conditions
				(true = fixed, false = free)
				1: x-axis minima boundary condition
				2: x-axis maxima boundary condition
				3: y-axis minima boundary condition
				4: y-axis maxima boundary condition
		output:
			X_m = {
				sin((m + 1)xπ / X),		dirichlet boundary condition
				cos(mxπ / X),			neumann boundary condition
				sin((m + 0.5)xπ / X),	mixed boundary conditions
				| m ∈ [0, ∞)
			}
			Y_n = {
				sin((n + 1)yπ / Y),		dirichlet boundary condition
				cos(nyπ / Y),			neumann boundary condition
				sin((n + 0.5)yπ / Y),	mixed boundary conditions
				| n ∈ [0, ∞)
			}
			U_xy = { X_m * Y_n | U ∈ ℝ^2 }
		*/

		T::Matrix_2D U(X, T::Matrix_1D(Y, 0.));
		const double pi_x = std::numbers::pi / (X - 1);
		const double pi_y = std::numbers::pi / (Y - 1);
		// boundary condition lambda
		auto BCLambda = [](const double mode,
						   const double scalar,
						   const bool minima,
						   const bool maxima) -> double {
			if (minima && maxima) {
				// dirichlet boundary
				return std::sin((mode + 1.) * scalar);
			} else if (!minima && !maxima) {
				// neumann boundary
				return std::cos(mode * scalar);
			} else {
				// mixed boundary
				return std::sin((mode + 0.5) * scalar);
			}
		};
		// produce domain U
		for (std::size_t x = 0; x < X; x++) {
			double X_m = BCLambda(m, x * pi_x, boundary_conditions[0], boundary_conditions[1]);
			for (std::size_t y = 0; y < Y; y++) {
				U[x][y] =
					X_m * BCLambda(n, y * pi_y, boundary_conditions[2], boundary_conditions[3]);
			}
		}
		return U;
	}

	inline T::Matrix_2D rectangularSeries(
		const std::size_t& M,
		const std::size_t& N,
		const double& epsilon,
		const std::array<bool, 4> boundary_conditions = {true, true, true, true}
	) {
		/*
		Calculate the eigenvalues of a 2-dimensional rectangular domain.
		input:
			M = number of modes across the Mth axis
			N = number of modes across the Nth axis
			epsilon = aspect ratio of the rectangle
			boundary_conditions = boolean array indicating the boundary conditions
				(true = fixed, false = free)
				1: x-axis minima boundary condition
				2: x-axis maxima boundary condition
				3: y-axis minima boundary condition
				4: y-axis maxima boundary condition
		output:
			X_m = {
				(m + 1)^2 / Є,			dirichlet boundary condition
				m^2 / Є,				neumann boundary condition
				(m + 0.5)^2 / Є,		mixed boundary conditions
				| m ∈ [0, M)
			}
			Y_n = {
				(n + 1)^2 * Є,			dirichlet boundary condition
				n^2 * Є,				neumann boundary condition
				(n + 0.5)^2 * Є,		mixed boundary conditions
				| n ∈ [0, N)
			}
			λ_mn = { √(X_m + Y_n) | λ ∈ ℝ }
		*/

		T::Matrix_2D S(M, T::Matrix_1D(N, 0.));
		const double epsilon_recip = 1. / epsilon;
		// boundary condition lambda
		auto BCLambda = [](const std::size_t i,
						   const double scalar,
						   const bool minima,
						   const bool maxima) -> double {
			if (minima && maxima) {
				// dirichlet boundary
				return (i + 1.) * (i + 1.) * scalar;
			} else if (!minima && !maxima) {
				// neumann boundary
				return i * i * scalar;
			} else {
				// mixed boundary
				return (i + 0.5) * (i + 0.5) * scalar;
			}
		};
		// produce the eigenvalues
		for (std::size_t m = 0; m < M; m++) {
			double m_hat =
				BCLambda(m, epsilon_recip, boundary_conditions[0], boundary_conditions[1]);
			for (std::size_t n = 0; n < N; n++) {
				S[m][n] = std::sqrt(
					m_hat + BCLambda(n, epsilon, boundary_conditions[2], boundary_conditions[3])
				);
			}
		}
		return S;
	}

}
