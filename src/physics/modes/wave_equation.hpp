/*
Functions relating to the wave equation.
*/

#pragma once

// core
#include <algorithm>	// max
#include <math.h>
#include <numbers>
#include <vector>
using namespace std::numbers;

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	inline T::Matrix_1D WaveEquationWaveform2D(
		T::Matrix_2D F,
		const T::Matrix_2D& A,
		const double& d,
		const double& k,
		const unsigned long& T
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
			waveform = W[t] ∈ A * e^dt * sin(ωt) / max(A) * NM
		*/

		T::Matrix_1D waveform(T);
		const unsigned long N = F.size();
		const unsigned long M = F[0].size();
		double A_max_NM = 0.;
		for (unsigned long n = 0; n < N; n++) {
			for (unsigned long m = 0; m < M; m++) {
				// calculate A_max and transform F into ω
				A_max_NM = std::max(A_max_NM, A[n][m]);
				F[n][m] *= 2 * pi * k;
			}
		}
		A_max_NM *= N * M;
		for (unsigned long t = 0; t < T; t++) {
			double d_t = pow(e, t * d);
			for (unsigned long n = 0; n < N; n++) {
				for (unsigned long m = 0; m < M; m++) {
					// 2009 - Bilbao, pp.65-66
					// 2016 - Chaigne & Kergomard, p.154
					waveform[t] += A[n][m] * d_t * sin(t * F[n][m]) / A_max_NM;
				}
			}
		}
		return waveform;
	}

}
