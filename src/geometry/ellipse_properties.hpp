/*
Utility functions for working with polygons.
*/

#pragma once

// core
#include <math.h>
#include <numbers>
#include <vector>
using namespace std::numbers;

// src
#include "../types.hpp"
namespace T = kac_core::types;

namespace kac_core::geometry {

	inline double circleArea(const double& r) {
		/*
		Archimedes area for a circle.
		input:
			r = radius
		*/
		return pi * pow(r, 2.0);
	}

}