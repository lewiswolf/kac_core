/*
Functions for calculating an audio waveform using a finite difference time
domain method.
*/

#pragma once

// core
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

// src
#include "../../types.hpp"
namespace T = kac_core::types;

namespace kac_core::physics {

	inline T::Matrix_1D FDTDWaveform1D(
		T::Matrix_1D u_0,
		T::Matrix_1D u_1,
		const double& c_0,
		const double& c_1,
		const double& c_2,
		const std::size_t& T,
		const double& w
	) {
		/*
		Generates a waveform using a 1-dimensional FDTD scheme.
		input:
			u_0 = initial fdtd grid at t = 0.
			u_1 = initial fdtd grid at t = 1.
			c_0 = first fdtd coefficient related to the decay term and the
				courant number.
			c_1 = second fdtd coefficient related to the decay term and the
				courant number.
			c_2 = third fdtd coefficient related to the decay term.
			T = length of simulation in samples.
			w = the coordinate at which the waveform is sampled ∈ ℝ^1, [0. 1.].
		output:
			waveform = W[n + 1] ∈ (λ ** 2)(
				u_n_x+1_y + u_n_x-1_y + u_n_x_y+1 + u_n_x_y-1
			) + 2(1 - 2(λ ** 2))u_n_x_y - d(u_n-1_x_y) ∀ u ∈ R^2
		*/

		// handle errors
		if (u_0.size() != u_1.size()) {
			throw std::invalid_argument("u_0 and u_1 differ in size.");
		}
		// lambda for sampling the 2D matrix using bilinear interpolation.
		const std::size_t X = u_0.size() - 2;
		const std::size_t x_0 = std::floor(w * X);
		const double a = w * X - x_0;
		const double coef_0 = (1. - a);
		const auto linearInterpolation = [=](const T::Matrix_1D& u) {
			return coef_0 * u[x_0] + a * u[x_0 + 1];
		};
		// initialise output
		T::Matrix_1D waveform(T);
		waveform[0] = linearInterpolation(u_0);
		waveform[1] = linearInterpolation(u_1);
		// lambda for the update equation
		const auto FDTDUpdate1D = [=](T::Matrix_1D& u_a, const T::Matrix_1D& u_b) {
			for (std::size_t x = 1; x <= X; x++) {
				// update in place
				u_a[x] = c_0 * (u_b[x + 1] + u_b[x - 1]) + c_1 * u_b[x] - c_2 * u_a[x];
			};
		};
		// main loop
		for (std::size_t t = 2; t < T; t++) {
			// maintain memory efficiency - only two matrices need to be in memory at any time
			FDTDUpdate1D(u_0, u_1);
			waveform[t] = linearInterpolation(u_0);
			std::swap(u_0, u_1);
		}
		// normalise
		double max_x = 0.;
		for (double& x: waveform) max_x = std::max(max_x, std::abs(x));
		if (max_x != 0.) {
			for (double& x: waveform) x /= max_x;
		}
		return waveform;
	}

	inline T::Matrix_1D FDTDWaveform2D(
		T::Matrix_2D u_0,
		T::Matrix_2D u_1,
		const T::BooleanImage_2D& B,
		const double& c_0,
		const double& c_1,
		const double& c_2,
		const std::size_t& T,
		const T::Point& w
	) {
		/*
		Generates a waveform using a 2 dimensional FDTD scheme.
		input:
			u_0 = initial fdtd grid at t = 0.
			u_1 = initial fdtd grid at t = 1.
			B = boundary conditions.
			c_0 = first fdtd coefficient related to the decay term and the
				courant number.
			c_1 = second fdtd coefficient related to the decay term and the
				courant number.
			c_2 = third fdtd coefficient related to the decay term.
			T = length of simulation in samples.
			w = the coordinate at which the waveform is sampled ∈ ℝ^2, [0. 1.].
		output:
			waveform = W[n + 1] ∈ (λ ** 2)(
				u_n_x+1_y + u_n_x_y+1 + u_n_x-1_y + u_n_x_y-1
			) + 2(1 - 2(λ ** 2))u_n_x_y - d(u_n-1_x_y) ∀ u ∈ R^2
		*/

		// handle errors
		if (u_0.size() != u_1.size() || u_0[0].size() != u_1[0].size()) {
			throw std::invalid_argument("u_0 and u_1 differ in size.");
		}
		if (u_0.size() != B.size() || u_0[0].size() != B[0].size()) {
			throw std::invalid_argument("u_0 and B differ in size.");
		}
		// lambda for sampling the 2D matrix using bilinear interpolation.
		const std::size_t x_0 = std::floor(w.x * (u_0.size() - 2));
		const std::size_t y_0 = std::floor(w.y * (u_0[0].size() - 2));
		const double a = w.x * (u_0.size() - 2) - x_0;
		const double b = w.y * (u_0[0].size() - 2) - y_0;
		const double coef_0 = (1. - a) * (1 - b);
		const double coef_1 = (1. - a) * b;
		const double coef_2 = a * (1. - b);
		const double coef_3 = a * b;
		const auto bilinearInterpolation = [=](const T::Matrix_2D& u) {
			return coef_0 * u[x_0][y_0] + coef_1 * u[x_0][y_0 + 1] + coef_2 * u[x_0 + 1][y_0]
				 + coef_3 * u[x_0 + 1][y_0 + 1];
		};
		// initialise output
		T::Matrix_1D waveform(T);
		waveform[0] = bilinearInterpolation(u_0);
		waveform[1] = bilinearInterpolation(u_1);
		// for efficiency, calculate the loop range relative to dirichlet boundary conditions
		const std::size_t B_dim_X = B.size();
		const std::size_t B_dim_Y = B[0].size();
		std::array<std::size_t, 2> x_range = {B_dim_X, 0};
		std::array<std::size_t, 2> y_range = {B_dim_Y, 0};
		// forward loop to find the first ones
		for (std::size_t x = 1; x < B_dim_X - 1; x++) {
			for (std::size_t y = 1; y < B_dim_Y - 1; y++) {
				if (B[x][y]) {
					x_range[0] = std::min(x, x_range[0]);
					y_range[0] = std::min(y, y_range[0]);
					break;
				}
			}
		}
		// backwards loop to find the last ones
		for (std::size_t x = B_dim_X - 2; x > 0; x--) {
			for (std::size_t y = B_dim_Y - 2; y > 0; y--) {
				if (B[x][y]) {
					x_range[1] = std::max(x, x_range[1]);
					y_range[1] = std::max(y, y_range[1]);
					break;
				}
			}
		}
		// lambda for the update equation
		const auto FDTDUpdate2D = [=](T::Matrix_2D& u_a, const T::Matrix_2D& u_b) {
			for (std::size_t x = x_range[0]; x <= x_range[1]; x++) {
				for (std::size_t y = y_range[0]; y <= y_range[1]; y++) {
					// dirichlet boundary conditions
					if (B[x][y]) {
						// update in place
						u_a[x][y] =
							c_0 * (u_b[x + 1][y] + u_b[x][y + 1] + u_b[x - 1][y] + u_b[x][y - 1])
							+ c_1 * u_b[x][y] - c_2 * u_a[x][y];
					}
				};
			}
		};
		// main loop
		for (std::size_t t = 2; t < T; t++) {
			// maintain memory efficiency - only two matrices need to be in memory at any time
			FDTDUpdate2D(u_0, u_1);
			waveform[t] = bilinearInterpolation(u_0);
			std::swap(u_0, u_1);
		}
		// normalise
		double max_x = 0.;
		for (double& x: waveform) max_x = std::max(max_x, std::abs(x));
		if (max_x != 0.) {
			for (double& x: waveform) x /= max_x;
		}
		return waveform;
	}

	inline T::Matrix_1D FDTDUpdate1D(
		T::Matrix_1D& u_0,
		const T::Matrix_1D& u_1,
		const double& c_0,
		const double& c_1,
		const double& c_2
	) {
		/*
		2-dimensional FDTD update equation.
		input:
			u_0 = initial fdtd grid at t = -1.
			u_1 = initial fdtd grid at t = 0.
			c_0 = first fdtd coefficient related to the decay term and the courant number.
			c_1 = second fdtd coefficient related to the decay term and the courant number.
			c_2 = third fdtd coefficient related to the decay term.
		output:
			u = c_0 * (u_1_x+1 + u_1_x-1) + c_1 * u_1_x - c_2 * u_0_x
		*/

		const std::size_t X = u_0.size() - 1;
		for (std::size_t x = 1; x < X; x++) {
			// update in place
			u_0[x] = c_0 * (u_1[x + 1] + u_1[x - 1]) + c_1 * u_1[x] - c_2 * u_0[x];
		}
		return u_0;
	}

	inline T::Matrix_2D FDTDUpdate2D(
		T::Matrix_2D& u_0,
		const T::Matrix_2D& u_1,
		const T::BooleanImage_2D& B,
		const double& c_0,
		const double& c_1,
		const double& c_2,
		const std::array<std::size_t, 2>& x_range,
		const std::array<std::size_t, 2>& y_range
	) {
		/*
		2-dimensional FDTD update equation.
		input:
			u_0 = initial fdtd grid at t = -1.
			u_1 = initial fdtd grid at t = 0.
			B = boundary conditions.
			c_0 = first fdtd coefficient related to the decay term and the courant number.
			c_1 = second fdtd coefficient related to the decay term and the courant number.
			c_2 = third fdtd coefficient related to the decay term.
			x_range = range across the x-axis of the boundary condition (for optimisation).
			y_range = range across the y-axis of the boundary condition (for optimisation).
		output:
			u = c_0 * (u_1_x+1_y + u_1_x_y+1 + u_1_x-1_y + u_1_x_y-1)
				+ c_1 * u_1_x_y - c_2 * u_0_x_y
		*/

		for (std::size_t x = x_range[0]; x <= x_range[1]; x++) {
			for (std::size_t y = y_range[0]; y <= y_range[1]; y++) {
				// dirichlet boundary conditions
				if (B[x][y]) {
					// update in place
					u_0[x][y] =
						c_0 * (u_1[x + 1][y] + u_1[x][y + 1] + u_1[x - 1][y] + u_1[x][y - 1])
						+ c_1 * u_1[x][y] - c_2 * u_0[x][y];
				}
			};
		}
		return u_0;
	}

}
