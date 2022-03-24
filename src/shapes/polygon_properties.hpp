bool isConvex(Vertices const& v) {
	/*
	Tests whether or not a given array of vertices forms a convex polygon. This
	is achieved using the resultant sign of the cross product for each vertex:
		[(x_i - x_i-1), (y_i - y_i-1)] x [(x_i+1 - x_i), (y_i+1 - y_i)]
	See => http://paulbourke.net/geometry/polygonmesh/ 'Determining whether or
	not a polygon (2D) has its vertices ordered clockwise or counter-clockwise'.
	*/

	// anonymous function declaration of cross product - z component early, see
	// np.cross =>
	// https://numpy.org/doc/stable/reference/generated/numpy.cross.html
	std::function<double(Point, Point, Point)> crossProductZ =
		[&](Point p, Point p_plus, Point p_minus) {
			return (p.x - p_minus.x) * (p_plus.y - p.y) -
				   (p_plus.x - p.x) * (p.y - p_minus.y);
		};

	// determine the direction of the initial point using the cross product.
	bool clockwise = crossProductZ(v[0], v[1], v[v.size() - 1]) < 0;

	// loop over remaining points
	for (unsigned int i = 1; i < v.size(); i++) {
		if (crossProductZ(v[i], v[(i + 1) % v.size()], v[i - 1]) < 0 !=
			clockwise) {
			return false;
		}
	}
	return true;
}