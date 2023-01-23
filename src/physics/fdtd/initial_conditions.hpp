/*
Functions for producing a raised cosine transform for different dimensionalities.
*/

#pragma once

// core
#include <math.h>
#include <numbers>
#include <vector>
using namespace std::numbers;

// src
#include "../../types.hpp"
using namespace kac_core::types;

namespace kac_core::physics {

	Matrix_1D raisedCosine1D(const int& size, const double& mu, const double& sigma) {
		/*
		Calculate a one dimensional raised cosine distribution.
		See Bilbao, S. (2009) Numerical Sound Synthesis p.121.
		input:
			size = the size of the matrix.
			μ = a cartesian point representing the maxima of the cosine.
			σ = variance.
		output:
			{
				(1 + cos(π(x - μ) / σ)) / 2,	|x - μ| ≤ σ
				0,								|x - μ| > σ
			}
		*/

		Matrix_1D raised_cosine(size);
		for (int x = 0; x < size; x++) {
			double x_diff = fabs(x - mu);
			if (x_diff <= sigma) {
				raised_cosine[x] = 0.5 * (1 + cos(pi * x_diff / sigma));
			}
		}
		return raised_cosine;
	}

	Matrix_2D raisedCosine2D(
		const int& size_X,
		const int& size_Y,
		const double& mu_x,
		const double& mu_y,
		const double& sigma
	) {
		/*
		Calculate a two dimensional raised cosine distribution.
		See Bilbao, S. (2009) Numerical Sound Synthesis p.306.
		input:
			size = the size of the matrix.
			μ = a cartesian point representing the maxima of the cosine.
			σ = variance.
		output:
			l2_norm = ((x - mu_x)^2 + (y - mu_y)^2)^0.5
			{
				(1 + cos(π(l2_norm) / σ)) / 2,	|l2_norm| ≤ σ
				0,								|l2_norm| > σ
			}
		*/

		Matrix_2D raised_cosine(size_X, Matrix_1D(size_Y, 0));
		for (int x = 0; x < size_X; x++) {
			for (int y = 0; y < size_Y; y++) {
				double l2_norm = sqrt(pow((x - mu_x), 2) + pow((y - mu_y), 2));
				if (l2_norm <= sigma) {
					raised_cosine[x][y] = 0.5 * (1 + cos(pi * l2_norm / sigma));
				}
			}
		}
		return raised_cosine;
	}

	Matrix_1D
	raisedTriangle1D(const int& size, const double& mu, const double& a, const double& b) {
		/*
		Calculate a one dimensional triangular distribution.
		See Bilbao, S. (2009) Numerical Sound Synthesis p.121.
		input:
			size = the size of the matrix.
			μ = a cartesian point representing the maxima of the triangle.
			a = minimum x value for the distribution.
			b = maximum x value for the distribution.
		output:
			Λ(x) = {
				0,								x < a
				(x - a) / (μ - a),				a ≤ x ≤ μ
				1. - (x - μ) / (b - μ),			μ < x ≤ b
				0,								x > a
			}
		*/

		Matrix_1D triangle(size);
		for (int x = 0; x < size; x++) {
			triangle[x] = a <= x && x <= mu ? double(x - a) / double(mu - a)
						: mu < x && x <= b	? 1.0 - double(x - mu) / double(b - mu)
											: 0.0;
		}
		return triangle;
	}

	Matrix_2D raisedTriangle2D(
		const int& size_X,
		const int& size_Y,
		const double& mu_x,
		const double& mu_y,
		const double& x_a,
		const double& x_b,
		const double& y_a,
		const double& y_b
	) {
		/*
		Calculate a two dimensional triangular distribution.
		See https://reference.wolfram.com/language/ref/UnitTriangle.html
		input:
			size = the size of the matrix.
			μ = a cartesian point representing the maxima of the triangle.
			x_ab = minimum and maximum x value for the distribution.
			y_ab = minimum and maximum y value for the distribution.
		output:
			Λ(x, y) = Λ(x) * Λ(y)
		*/

		Matrix_1D y_t = raisedTriangle1D(size_Y, mu_y, y_a, y_b);
		Matrix_2D triangle(size_X, Matrix_1D(size_Y, 0));
		for (int x = 0; x < size_X; x++) {
			double x_t = x_a <= x && x <= mu_x ? double(x - x_a) / double(mu_x - x_a)
					   : mu_x < x && x <= x_b  ? 1.0 - double(x - mu_x) / double(x_b - mu_x)
											   : 0.0;
			for (int y = 0; y < size_Y; y++) { triangle[x][y] = x_t * y_t[y]; }
		}
		return triangle;
	}

}