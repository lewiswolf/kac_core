/*
Functions for calculating the linear approximation of a 1-Dimensional wave
equation.
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

	inline T::Matrix_1D linearAmplitudes(const double& x, const unsigned long& N) {
		/*
		Calculate the amplitudes of the 1D eigenmodes relative to a strike location.
		input:
			x = strike location
			N = number of modes
		output:
			A = { abs(sin(nxπ)) | a ∈ ℝ, 0 < n <= N }
		*/

		T::Matrix_1D A(N);
		double x_pi = x * pi;
		for (unsigned long n = 0; n < N; n++) { A[n] = abs(sin((n + 1) * x_pi)); };
		return A;
	}

	inline T::Matrix_1D linearCymatics(
		const double& n, const unsigned long& H, const std::pair<bool, bool>& boundary_conditions
	) {
		/*
		Produce the 1D continuos cymatic diagram for a particular mode of a linear domain.
		input:
			n = nth modal index
			H = length of the X axis
			boundary_conditions = boolean pair indicating the boundary conditions
				- first = left boundary condition (true = fixed, false = free)
				- second = right boundary condition (true = fixed, false = free)
		output:
			M = {
				sin(nπx/H),			(true, true)
				cos(nπx/H),			(false, false)
				sin((n-0.5)πx/H),	(true, false) || (false, true)
			}
		*/

		T::Matrix_1D M(H);
		double omega = 0.0;
		if (boundary_conditions.first) {
			if (boundary_conditions.second) {
				// fixed left and right boundary condition
				omega = n * pi / H;
				for (unsigned long x = 0; x < H; x++) { M[x] = sin(omega * x); }
			} else {
				// fixed left, free right boundary condition
				omega = (n - 0.5) * pi / H;
				for (unsigned long x = 0; x < H; x++) { M[x] = sin(omega * x); }
			}
		} else {
			if (boundary_conditions.second) {
				// fixed right, free left boundary condition
				omega = (n - 0.5) * pi / H;
				for (unsigned long x = 0; x < H; x++) { M[x] = sin(omega * (H - 1 - x)); }
			} else {
				// free left and right boundary condition
				omega = n * pi / H;
				for (unsigned long x = 0; x < H; x++) { M[x] = cos(omega * x); }
			}
		}
		return M;
	}

	inline T::Matrix_1D linearSeries(const unsigned long& N) {
		/*
		Calculate the the harmonic series.
		input:
			N = number of modes
		output:
			S = { n | s ∈ ℕ, 0 < n <= N }
		*/

		T::Matrix_1D S(N);
		for (unsigned long n = 0; n < N; n++) { S[n] = n + 1; };
		return S;
	}

}
