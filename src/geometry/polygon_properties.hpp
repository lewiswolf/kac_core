/*
Utility functions for working with polygons.
*/

#pragma once

// src
#include "../types.hpp"
using namespace kac_core::types;

namespace kac_core { namespace geometry {

	bool isColinear(const Point& a, const Point& b, const Point& c) {
		/*
		Determines whether or not a given set of three vertices are
		colinear.
		*/

		return (c.y - b.y) * (b.x - a.x) == (b.y - a.y) * (c.x - b.x);
	}

	bool isConvex(const Vertices& v) {
		/*
		Tests whether or not a given array of vertices forms a convex
		polygon. This is achieved using the resultant sign of the cross
		product for each vertex:
			[(x_i - x_i-1), (y_i - y_i-1)] x [(x_i+1 - x_i), (y_i+1 - y_i)]
		See => http://paulbourke.net/geometry/polygonmesh/ 'Determining
		whether or not a polygon (2D) has its vertices ordered clockwise or
		counter-clockwise'.
		*/

		// cross product - z component only, see np.cross =>
		// https://numpy.org/doc/stable/reference/generated/numpy.cross.html
		auto crossProductZ = [](Point p, Point p_plus, Point p_minus) {
			return (p.x - p_minus.x) * (p_plus.y - p.y)
				- (p_plus.x - p.x) * (p.y - p_minus.y);
		};
		// determine the direction of the initial point using the cross
		// product
		bool clockwise = crossProductZ(v[0], v[1], v[v.size() - 1]) < 0;
		// loop over remaining points
		for (unsigned int i = 1; i < v.size(); i++) {
			if (crossProductZ(v[i], v[(i + 1) % v.size()], v[i - 1]) < 0
				!= clockwise) {
				return false;
			}
		}
		return true;
	}

}}