/*
Functions for calculating the bessel functions and their zero crossings.
*/

#pragma once

// dependencies
#include <boost/math/special_functions/bessel.hpp>

namespace geometry {

	double besselJ(const int& n, const double& x) {
		/*
		Calculates the bessel function of the first kind J_n(x).
		Adapted from void bess()
		http://www.falstad.com/circosc-java/CircOsc.java
		input:
			n = order of the bessel function
			x = x coordinate
		output:
			y = J_n(x) | y ∈ ℝ
		*/
		return boost::math::cyl_bessel_j(n, x);
	}

	double besselJZero(const double& n, const int& m) {
		/*
		Calculates the mth zero crossing of the bessel functions of the first
		kind. Adapted from `double zeroj()`, see =>
		http://www.falstad.com/circosc-java/CircOsc.java
		input:
			n = order of the bessel function
			m = mth zero
		output:
			z_nm = mth zero crossing of J_n() | z_mn ∈ ℝ
		*/
		return boost::math::cyl_bessel_j_zero(n, m);
	}

}