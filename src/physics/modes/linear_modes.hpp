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

	inline T::Matrix_1D linearCymatics(const double& n, const unsigned long& H) {
		/*
		Produce the 1D continuos cymatic diagram for a particular mode of a linear domain.
		input:
			n = nth modal index
			H = length of the X axis
		output:
			M = cos(nπx/H) ≈ 0
		*/

		T::Matrix_1D M(H);
		for (unsigned long x = 0; x < H; x++) { M[x] = cos(n * pi * x / H); }
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
