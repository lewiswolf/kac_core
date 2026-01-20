/*
Functions for calculating the linear approximation of a 1-dimensional wave equation.
*/

#pragma once

// core
#include <cmath>
#include <numbers>
#include <vector>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	inline T::Matrix_1D linearAmplitudes(
		const double& x,
		const std::size_t& N,
		const std::array<bool, 2> boundary_conditions = {true, true}
	) {
		/*
		Calculate the spatial eigenfunction of a 1-dimensional domain relative to a strike location.
		input:
			x = strike location
			N = number of modes
			boundary_conditions = boolean array indicating the boundary conditions
				(true = fixed, false = free)
				1: x-axis minima boundary condition
				2: x-axis maxima boundary condition
		output:
			α_n = {
				sin((n + 1)πx),			dirichlet boundary condition
				cos(nπx),				neumann boundary condition
				sin((n + 0.5)πx),		mixed boundary conditions
				| α ∈ ℝ, n ∈ [0, N)
			}
		*/

		T::Matrix_1D A(N, 0.);
		const double x_pi = x * std::numbers::pi;
		if (boundary_conditions[0] && boundary_conditions[0]) {
			// dirichlet boundary
			for (std::size_t n = 0; n < N; n++) { A[n] = std::sin((n + 1.) * x_pi); }
		} else if (!boundary_conditions[0] && !boundary_conditions[0]) {
			// neumann boundary
			for (std::size_t n = 0; n < N; n++) { A[n] = std::cos(n * x_pi); }
		} else {
			// mixed boundary
			for (std::size_t n = 0; n < N; n++) { A[n] = std::sin((n + 0.5) * x_pi); }
		};
		return A;
	}

	inline T::Matrix_1D linearCymatics(
		const double& n,
		const std::size_t& X,
		const std::array<bool, 2> boundary_conditions = {true, true}
	) {
		/*
		Produce the cymatic diagram of a 1-dimensional domain for a particular mode λ_n.
		input:
			n = nth modal index
			X = length of the X axis
			boundary_conditions = boolean array indicating the boundary conditions
				(true = fixed, false = free)
				1: x-axis minima boundary condition
				2: x-axis maxima boundary condition
		output:
			U_x = {
				sin((n + 1) πx/H),		dirichlet boundary condition
				cos(nπx/H),				neumann boundary condition
				sin((n + 0.5)πx/H),		mixed boundary conditions
				| U ∈ ℝ^1, n ∈ [0, ∞)
			}
		*/

		T::Matrix_1D U(X, 0.);
		double omega = std::numbers::pi / X;
		if (boundary_conditions[0] && boundary_conditions[1]) {
			// dirichlet boundary
			omega *= (n + 1.);
			for (std::size_t x = 0; x < X; x++) { U[x] = std::sin(omega * x); }
		} else if (!boundary_conditions[0] && !boundary_conditions[1]) {
			// neumann boundary
			omega *= n;
			for (std::size_t x = 0; x < X; x++) { U[x] = std::cos(omega * x); }
		} else {
			// mixed boundary
			omega *= (n + 0.5);
			for (std::size_t x = 0; x < X; x++) { U[x] = std::sin(omega * x); }
		}
		return U;
	}

	inline T::Matrix_1D linearSeries(
		const std::size_t& N, const std::array<bool, 2> boundary_conditions = {true, true}
	) {
		/*
		Calculate the eigenvalues of a 1-dimensional domain.
		input:
			N = number of modes
			boundary_conditions = boolean array indicating the boundary conditions
				(true = fixed, false = free)
				1: x-axis minima boundary condition
				2: x-axis maxima boundary condition
		output:
			λ_n = {
				n + 1,					dirichlet boundary condition
				n,						neumann boundary condition
				n + 0.5,				mixed boundary conditions
				| λ ∈ ℝ, n ∈ [0, N)
			}
		*/

		T::Matrix_1D S(N, 0.);
		if (boundary_conditions[0] && boundary_conditions[1]) {
			// dirichlet boundary
			for (std::size_t n = 0; n < N; n++) { S[n] = n + 1.; }
		} else if (!boundary_conditions[0] && !boundary_conditions[1]) {
			// neumann boundary
			for (std::size_t n = 0; n < N; n++) { S[n] = n; }
		} else {
			// mixed boundary
			for (std::size_t n = 0; n < N; n++) { S[n] = n + 0.5; }
		};
		return S;
	}

}
