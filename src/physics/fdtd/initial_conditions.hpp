/*
Functions for producing a raised cosine transform for different dimensionalities.
*/

#pragma once

// core
#include <math.h>
#include <numbers>
#include <stdexcept>
#include <vector>
using namespace std::numbers;

// src
#include "../../types.hpp"
using namespace kac_core::types;

namespace kac_core::physics {

	inline Matrix_1D
	raisedCosine1D(const double& mu, const double& sigma, const std::size_t& size) {
		/*
		Calculate a one dimensional raised cosine distribution, normalised to a unit interval.
		See Bilbao, S. (2009) Numerical Sound Synthesis p.121.
		input:
			μ = a normalised point representing the maxima of the distribution ∈ [0, 1].
			σ = normalised variance ∈ (0, ∞].
			size = the size of the matrix.
		output:
			{
				(1 + cos(π(x - μ) / σ)) / 2,	|x - μ| ≤ σ
				0,								|x - μ| > σ
			}
		*/

		Matrix_1D raised_cosine(size, 0.);
		if (sigma > 0.) {
			const double inv_X = (size > 1) ? 1. / static_cast<double>(size - 1) : 0.;
			for (std::size_t i = 0; i < size; i++) {
				double x_diff = abs((static_cast<double>(i) * inv_X) - mu);
				if (x_diff <= sigma) {
					raised_cosine[i] = 0.5 * (1. + cos(pi * x_diff / sigma));
				}
			}
		}
		return raised_cosine;
	}

	inline Matrix_2D raisedCosine2D(
		const T::Point& mu,
		const double& sigma,
		const std::size_t& size_X,
		const std::size_t& size_Y
	) {
		/*
		Calculate a two dimensional raised cosine distribution, normalised to a unit interval.
		See Bilbao, S. (2009) Numerical Sound Synthesis p.306.
		input:
			μ = a normalised point representing the maxima of the distribution ∈ [0, 1].
			σ = normalised variance ∈ (0, ∞].
			size = the size of the matrix.
		output:
			l2_norm = ((x - mu_x)^2 + (y - mu_y)^2)^0.5
			{
				(1 + cos(π(l2_norm) / σ)) / 2,	|l2_norm| ≤ σ
				0,								|l2_norm| > σ
			}
		*/

		Matrix_2D raised_cosine(size_X, Matrix_1D(size_Y, 0.));
		if (sigma > 0.) {
			const double inv_X = (size_X > 1) ? 1. / static_cast<double>(size_X - 1) : 0.;
			const double inv_Y = (size_Y > 1) ? 1. / static_cast<double>(size_Y - 1) : 0.;
			for (std::size_t i = 0; i < size_X; i++) {
				double x = (static_cast<double>(i) * inv_X);
				for (std::size_t j = 0; j < size_Y; j++) {
					double l2_norm = std::hypot(x - mu.x, (static_cast<double>(j) * inv_Y) - mu.y);
					if (l2_norm <= sigma) {
						raised_cosine[i][j] = 0.5 * (1. + cos(pi * l2_norm / sigma));
					}
				}
			}
		}
		return raised_cosine;
	}

	inline Matrix_1D raisedTriangle1D(
		const double& mu, const double& x_a, const double& x_b, const std::size_t& size
	) {
		/*
		Calculate a one dimensional triangular distribution.
		See Bilbao, S. (2009) Numerical Sound Synthesis p.121.
		input:
			μ = a normalised point representing the maxima of the distribution ∈ [0, 1].
			x_a = segment length of distribution such that a = μ - x_a.
			x_b = segment length of distribution such that b = μ - x_b.
			size = the size of the matrix.
		output:
			Λ(x) = {
				0,								x < a
				(x - a) / (μ - a),				a ≤ x ≤ μ
				1. - (x - μ) / (b - μ),			μ < x ≤ b
				0,								x > b
			}
		*/

		Matrix_1D triangle(size, 0.);
		const double a = mu - std::max(x_a, 0.);
		const double b = mu + std::max(x_b, 0.);
		const double inv_X = (size > 1) ? 1. / static_cast<double>(size - 1) : 0.;
		for (std::size_t i = 0; i < size; i++) {
			double x = (static_cast<double>(i) * inv_X);
			triangle[i] = a <= x && x <= mu ? (x - a) / (mu - a)
						: mu < x && x <= b	? 1. - (x - mu) / (b - mu)
											: 0.;
		}
		return triangle;
	}

	inline Matrix_2D raisedTriangle2D(
		const T::Point& mu,
		const double& x_a,
		const double& x_b,
		const double& y_a,
		const double& y_b,
		const std::size_t& size_X,
		const std::size_t& size_Y
	) {
		/*
		Calculate a two dimensional triangular distribution.
		See https://reference.wolfram.com/language/ref/UnitTriangle.html
		input:
			μ = a normalised point representing the maxima of the distribution ∈ [0, 1].
			x_a = segment length of horizontal distribution such that a = μ - x_a.
			x_b = segment length of horizontal distribution such that b = μ - x_b.
			y_a = segment length of vertical distribution such that a = μ - y_a.
			y_b = segment length of vertical distribution such that b = μ - y_b.
			size = the size of the matrix.
		output:
			Λ(x, y) = Λ(x) * Λ(y)
		*/

		Matrix_2D triangle(size_X, Matrix_1D(size_Y, 0.));
		Matrix_1D x_t = raisedTriangle1D(mu.x, x_a, x_b, size_X);
		Matrix_1D y_t = raisedTriangle1D(mu.y, y_a, y_b, size_Y);
		for (std::size_t x = 0; x < size_X; x++) {
			for (std::size_t y = 0; y < size_Y; y++) { triangle[x][y] = x_t[x] * y_t[y]; }
		}
		return triangle;
	}

}
