/*
Functions relating to the wave equation.
*/

#pragma once

// core
#include <algorithm>	// max
#include <cmath>
#include <numbers>
#include <vector>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

T::Matrix_1D normalise(T::Matrix_1D& waveform) {
	double max_A = 0.;
	for (double& x: waveform) max_A = std::max(max_A, std::abs(x));
	if (max_A != 0.) {
		for (double& x: waveform) x /= max_A;
	}
	return waveform;
}

namespace kac_core::physics {

	inline T::Matrix_1D AdditiveSynthesis1D(
		const T::Matrix_1D& F,
		const T::Matrix_1D& alpha,
		const double& d,
		const double& k,
		const std::size_t& T
	) {
		/*
		Calculate a closed form solution to the 1D wave equation.
		input:
			F = frequencies (hertz)
			α = amplitudes ∈ [0, 1]
			d = decay
			k = sample length
			T = length of simulation
		output:
			waveform = W[t] ∈ e^dt * sin(ωt) * α
		*/

		T::Matrix_1D waveform(T);
		const std::size_t N = F.size();
		const double radians = 2. * std::numbers::pi * k;
		const double d_step = std::exp(d);
		double d_t = 1.;
		for (std::size_t t = 0; t < T; t++) {
			double t_radians = t * radians;
			for (std::size_t n = 0; n < N; n++) {
				// 2009 - Bilbao, pp.65-66
				// 2016 - Chaigne & Kergomard, p.154
				waveform[t] += std::sin(F[n] * t_radians) * alpha[n];
			}
			waveform[t] *= d_t;
			d_t *= d_step;
		}
		return normalise(waveform);
	}

	inline T::Matrix_1D AdditiveSynthesis2D(
		const T::Matrix_2D& F,
		const T::Matrix_2D& alpha,
		const double& d,
		const double& k,
		const std::size_t& T
	) {
		/*
		Calculate a closed form solution to the 2D wave equation.
		input:
			F = frequencies (hertz)
			α = amplitudes ∈ [0, 1]
			d = decay
			k = sample length
			T = length of simulation
		output:
			waveform = W[t] ∈ e^dt * sin(ωt) * α
		*/

		T::Matrix_1D waveform(T);
		const std::size_t N = F.size();
		const std::size_t M = F[0].size();
		const double radians = 2. * std::numbers::pi * k;
		const double d_step = std::exp(d);
		double d_t = 1.;
		for (std::size_t t = 0; t < T; t++) {
			double t_radians = t * radians;
			for (std::size_t n = 0; n < N; n++) {
				for (std::size_t m = 0; m < M; m++) {
					// 2009 - Bilbao, pp.65-66
					// 2016 - Chaigne & Kergomard, p.154
					waveform[t] += std::sin(F[n][m] * t_radians) * alpha[n][m];
				}
			}
			waveform[t] *= d_t;
			d_t *= d_step;
		}
		return normalise(waveform);
	}

	inline T::BooleanImage_1D
	ChladniPattern1D(const T::Matrix_1D& U, const double& tolerance = 0.1) {
		/*
		Produce a Chladni pattern from a 1-dimensional cymatic diagram.
		input:
			U = spatial eigenfunction ∈ [-1, 1]
			tolerance = thickness-dependent of the nodal lines
		output:
			B_x = abs(U_x) ≈ 0
		*/

		const std::size_t X = U.size();
		T::BooleanImage_1D B(X, 0);
		for (std::size_t x = 0; x < X; x++) { B[0][x] = std::abs(U[x]) < tolerance ? 1 : 0; }
		return B;
	}

	inline T::BooleanImage_2D
	ChladniPattern2D(const T::Matrix_2D& U, const double& tolerance = 0.1) {
		/*
		Produce a Chladni pattern from a 2-dimensional cymatic diagram.
		input:
			U = spatial eigenfunction ∈ [-1, 1]
			tolerance = thickness-dependent of the nodal lines
		output:
			B_xy = abs(U_xy) ≈ 0
		*/

		const std::size_t X = U.size();
		const std::size_t Y = U[0].size();
		T::BooleanImage_2D B(X, std::vector<short>(Y, 0));
		for (std::size_t x = 0; x < X; x++) {
			for (std::size_t y = 0; y < Y; y++) { B[x][y] = std::abs(U[x][y]) < tolerance ? 1 : 0; }
		}
		return B;
	}

}
