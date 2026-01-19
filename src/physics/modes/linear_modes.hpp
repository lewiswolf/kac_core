/*
Functions for calculating the linear approximation of a 1-Dimensional wave
equation.
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
		const std::pair<bool, bool> boundary_conditions = {true, true}
	) {
		/*
		Calculate the spatial eigenfunction of a 1D domain relative to a strike location.
		input:
			x = strike location
			N = number of modes
			boundary_conditions = boolean pair indicating the boundary conditions
				- first = left boundary condition (true = fixed, false = free)
				- second = right boundary condition (true = fixed, false = free)
		output:
			α_n = {
				sin((n + 1)πx),			dirichlet boundary conditions
				cos(nπx),				neumann boundary conditions
				sin((n + 0.5)πx),		mixed boundary conditions
				| α ∈ ℝ, 0 <= n < N
			}
		*/

		T::Matrix_1D A(N, 0.);
		const double x_pi = x * std::numbers::pi;
		if (boundary_conditions.first && boundary_conditions.second) {
			// dirichlet boundary
			for (std::size_t n = 0; n < N; n++) { A[n] = std::abs(std::sin((n + 1) * x_pi)); }
		} else if (!boundary_conditions.first && !boundary_conditions.second) {
			// neumann boundary
			for (std::size_t n = 0; n < N; n++) { A[n] = std::abs(std::cos(n * x_pi)); }
		} else {
			// mixed boundary
			for (std::size_t n = 0; n < N; n++) { A[n] = std::abs(std::sin((n + 0.5) * x_pi)); }
		};
		return A;
	}

	inline T::Matrix_1D linearCymatics(
		const double& n,
		const std::size_t& H,
		const std::pair<bool, bool> boundary_conditions = {true, true}
	) {
		/*
		Produce a cymatic diagram of a 1D domain for a particular mode n.
		input:
			n = nth modal index
			H = length of the X axis
			boundary_conditions = boolean pair indicating the boundary conditions
				- first = left boundary condition (true = fixed, false = free)
				- second = right boundary condition (true = fixed, false = free)
		output:
			U = {
				sin((n + 1) πx/H),		dirichlet boundary conditions
				cos(nπx/H),				neumann boundary conditions
				sin((n + 0.5)πx/H),		mixed boundary conditions
				| U ∈ ℝ^1, 0 <= n < N
			}
		*/

		T::Matrix_1D U(H, 0.);
		double omega = std::numbers::pi / H;
		if (boundary_conditions.first && boundary_conditions.second) {
			// dirichlet boundary
			omega *= (n + 1);
			for (std::size_t x = 0; x < H; x++) { U[x] = std::sin(omega * x); }
		} else if (!boundary_conditions.first && !boundary_conditions.second) {
			// neumann boundary
			omega *= n;
			for (std::size_t x = 0; x < H; x++) { U[x] = std::cos(omega * x); }
		} else {
			// mixed boundary
			omega *= (n + 0.5);
			for (std::size_t x = 0; x < H; x++) { U[x] = std::sin(omega * x); }
		}
		return U;
	}

	inline T::Matrix_1D linearSeries(
		const std::size_t& N, const std::pair<bool, bool> boundary_conditions = {true, true}
	) {
		/*
		Calculate the eigenvalues of a 1D domain.
		input:
			N = number of modes
			boundary_conditions = boolean pair indicating the boundary conditions
				- first = left boundary condition (true = fixed, false = free)
				- second = right boundary condition (true = fixed, false = free)
		output:
			λ_n = {
				n + 1,					dirichlet boundary conditions
				n,						neumann boundary conditions
				n + 0.5,				mixed boundary conditions
				| s ∈ ℝ, 0 <= n < N
			}
		*/

		T::Matrix_1D S(N, 0.);
		if (boundary_conditions.first && boundary_conditions.second) {
			// dirichlet boundary
			for (std::size_t n = 0; n < N; n++) { S[n] = n + 1; }
		} else if (!boundary_conditions.first && !boundary_conditions.second) {
			// neumann boundary
			for (std::size_t n = 0; n < N; n++) { S[n] = n; }
		} else {
			// mixed boundary
			for (std::size_t n = 0; n < N; n++) { S[n] = n + 0.5; }
		};
		return S;
	}

}
