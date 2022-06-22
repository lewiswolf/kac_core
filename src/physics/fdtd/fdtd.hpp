#pragma once

// core
#include <array>
#include <math.h>
#include <stdexcept>
#include <vector>

// src
#include "../../types.hpp"
using namespace kac_core::types;

namespace kac_core { namespace geometry {

	Matrix_1D FDTDWaveform2D(
		Matrix_2D u_0,
		Matrix_2D u_1,
		const std::vector<std::vector<char>>& B,
		const double& c_0,
		const double& c_1,
		const double& d,
		const int& T,
		Point w
	) {
		/*
		Generates a waveform using a 2 dimensional FDTD scheme.
		input:
			u_0 = initial fdtd grid at t = 0.
			u_1 = initial fdtd grid at t = 1.
			B = boundary conditions.
			c_0 = first fdtd coefficient related to the courant number.
			c_1 = second fdtd coefficient related to the courant number.
			d = decay coefficient.
			T = length of simulation in samples.
			w = the coordinate at which the waveform is sampled.
		output:
			waveform = W[n] ∈
				(λ ** 2)(
					u_n_x+1_y + u_n_x-1_y + u_n_x_y+1 + u_n_x_y-1
				) + 2(1 - 2(λ ** 2))u_n_x_y - d(u_n-1_x_y) ∀ u ∈ R^2
		*/

		// handle errors
		if (u_0.size() != u_1.size() || u_0[0].size() != u_1[0].size()) {
			throw std::invalid_argument("u_0 and u_1 differ in size.");
		}
		if (u_0.size() != B.size() || u_0[0].size() != B[0].size()) {
			throw std::invalid_argument("u_0 and B differ in size.");
		}
		// initialise output
		std::vector<double> waveform(T);
		waveform[0] = u_0[w.x][w.y];
		waveform[1] = u_1[w.x][w.y];
		// for efficiency, calculate the loop range relative to dirichlet
		// boundary conditions
		std::array<unsigned int, 2> x_range = {(unsigned int)B.size(), 0};
		std::array<unsigned int, 2> y_range = {(unsigned int)B[0].size(), 0};
		// forward loop to find the first ones
		for (unsigned int x = 0; x < B.size(); x++) {
			for (unsigned int y = 0; y < B[0].size(); y++) {
				if (B[x][y] == 1) {
					x_range[0] = x_range[0] > x ? x : x_range[0];
					y_range[0] = y_range[0] > y ? y : y_range[0];
					continue;
				}
			}
		}
		// backwards loop to find the last ones
		for (unsigned int x = B.size() - 1; x >= 0; x--) {
			for (unsigned int y = B[0].size() - 1; y >= 0; y--) {
				if (B[x][y] == 1) {
					x_range[1] = x_range[1] < x ? x : x_range[1];
					y_range[1] = y_range[1] < y ? y : y_range[1];
					continue;
				}
			}
		}

		// update lambda
		auto FDTDUpdate2D = [=](Matrix_2D& u_a, Matrix_2D& u_b) {
			for (unsigned int x = x_range[0]; x < x_range[1]; x++) {
				for (unsigned int y = y_range[0]; y < y_range[1]; y++) {
					// dirichlet boundary conditions
					if (B[x][y] != 0) {
						u_a[x][y] = (u_b[x][y + 1] + u_b[x + 1][y]
									 + u_b[x][y - 1] + u_b[x - 1][y])
								* c_0
							+ c_1 * u_b[x][y] - d * u_a[x][y];
					}
				};
			}
		};
		// main loop
		for (unsigned int t = 2; t < T; t++) {
			// branching maintains memory efficiency, meaning that only two
			// matrices need to be in memory at one time
			if ((t % 2) == 0) {
				FDTDUpdate2D(u_0, u_1);
				waveform[t] = u_0[w.x][w.y];
			} else {
				FDTDUpdate2D(u_1, u_0);
				waveform[t] = u_1[w.x][w.y];
			}
		}
		return waveform;
	}
}}