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

namespace kac_core::physics {

	inline T::Matrix_1D AdditiveSynthesis1D(
		const T::Matrix_1D& F,
		const T::Matrix_1D& alpha,
		const double& d,
		const double& k,
		const std::size_t& T
	) {
		/*
		Create a waveform of a 1-dimensional material using a physically informed
		representation of additive synthesis.
		input:
			F = frequencies (hertz)
			α = spatial eigenfunction ∈ [-1, 1]
			d = decay ∈ [0, ∞)
			k = sample length (ms)
			T = length of simulation (seconds)
		output:
			W[t] = ∑ e^dt * sin(f_n 2πkt) * α_n
		*/

		T::Matrix_1D waveform(T, 0.);
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
		// normalise
		double max_x = 0.;
		for (double& x: waveform) max_x = std::max(max_x, std::abs(x));
		if (max_x != 0.) {
			for (double& x: waveform) x /= max_x;
		}
		return waveform;
	}

	inline T::Matrix_1D AdditiveSynthesis2D(
		const T::Matrix_2D& F,
		const T::Matrix_2D& alpha,
		const double& d,
		const double& k,
		const std::size_t& T
	) {
		/*
		Create a waveform of a 2-dimensional material using a physically informed
		representation of additive synthesis.
		input:
			F = frequencies (hertz)
			α = amplitudes ∈ [-1, 1]
			d = decay ∈ [0, ∞)
			k = sample length (ms)
			T = length of simulation (seconds)
		output:
			W[t] = ∑ e^dt * sin(f_mn 2πkt) * α_mn
		*/

		T::Matrix_1D waveform(T, 0.);
		const std::size_t M = F.size();
		const std::size_t N = F[0].size();
		const double radians = 2. * std::numbers::pi * k;
		const double d_step = std::exp(d);
		double d_t = 1.;
		for (std::size_t t = 0; t < T; t++) {
			double t_radians = t * radians;
			for (std::size_t m = 0; m < M; m++) {
				for (std::size_t n = 0; n < N; n++) {
					// 2009 - Bilbao, pp.65-66
					// 2016 - Chaigne & Kergomard, p.154
					waveform[t] += std::sin(F[m][n] * t_radians) * alpha[m][n];
				}
			}
			waveform[t] *= d_t;
			d_t *= d_step;
		}
		// normalise
		double max_x = 0.;
		for (double& x: waveform) max_x = std::max(max_x, std::abs(x));
		if (max_x != 0.) {
			for (double& x: waveform) x /= max_x;
		}
		return waveform;
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
		for (std::size_t x = 0; x < X; x++) { B[x] = std::abs(U[x]) < tolerance ? 1 : 0; }
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
