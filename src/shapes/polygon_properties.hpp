bool isConvex(Vertices const& v) {
	/*
	Tests whether or not a given array of vertices forms a convex polygon. This
	is achieved using the resultant sign of the cross product for each vertex:
		[(x_i - x_i-1), (y_i - y_i-1)] x [(x_i+1 - x_i), (y_i+1 - y_i)]
	See => http://paulbourke.net/geometry/polygonmesh/ 'Determining whether or
	not a polygon (2D) has its vertices ordered clockwise or counter-clockwise'.
	*/

	// determine whether the initial point is clockwise using the above cross
	// product.
	bool clockwise = ((v[0].x - v[v.size() - 1].x) * (v[1].y - v[0].y) -
					  (v[1].x - v[0].x) * (v[0].y - v[v.size() - 1].y)) < 0;

	// loop over remaining points
	for (unsigned int i = 1; i < v.size(); i++) {
		if (((v[i].x - v[i - 1].x) * (v[(i + 1) % v.size()].y - v[i].y) -
			 (v[(i + 1) % v.size()].x - v[i].x) * (v[i].y - v[i - 1].y)) < 0 !=
			clockwise) {
			return false;
		}
	}
	return true;
}