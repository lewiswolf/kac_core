/*
Functions relating to the wave equation.
*/

#pragma once

// core
#define _USE_MATH_DEFINES
#include <algorithm>	// max
#include <math.h>
#include <vector>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	T::Matrix_1D WaveEquationWaveform2D(
		T::Matrix_2D F,
		const T::Matrix_2D& A,
		const double& d,
		const double& k,
		const unsigned int& T
	) {
		/*
		Calculate a closed form solution to the 2D wave equation.
		input:
			F = frequencies (hertz)
			A = amplitudes ∈ [0, 1]
			d = decay
			k = sample length
			T = length of simulation
		output:
			waveform = W[t] ∈ A * e^dt * sin(ωt) / NM * max(A)
		*/

		T::Matrix_1D waveform(T);
		unsigned int N = F.size();
		unsigned int M = F[0].size();
		double A_max = 0.0;
		for (unsigned int n = 0; n < N; n++) {
			for (unsigned int m = 0; m < M; m++) {
				// calculate A_max and transform F into ω
				A_max = std::max(A_max, A[n][m]);
				F[n][m] *= 2 * M_PI * k;
			}
		}
		for (unsigned int t = 0; t < T; t++) {
			double sum = 0.0;
			double d_t = pow(M_E, t * d);
			for (unsigned int n = 0; n < N; n++) {
				for (unsigned int m = 0; m < M; m++) {
					// 2009 - Bilbao, pp.65-66
					sum += A[n][m] * d_t * sin(t * F[n][m]) / (N * M * A_max);
				}
			}
			waveform[t] = sum;
		}
		return waveform;
	}

}