#include <algorithm>
#include <vector>

std::vector<double> calculateLinearSeries(int const& N) {
	/*
	Calculate the relationship between the modes of the harmonic series.

	params:
		N	-	number of modes
	*/

	std::vector<double> series(N);
	unsigned int i = 1;	 // cpp11 does not support declaring the iterator in the
						 // lambda function
	std::generate(series.begin(), series.end(), [&]() mutable {
		i++;
		return i;
	});
	return series;
}

std::vector<double> calculateLinearModes(double const& f0, int const& N) {
	/*
	Calculate the harmonic series of a given fundamental.

	params:
		f0	-	fundamental frequency
		N	-	number of modes
	*/

	std::vector<double> modes(N);
	unsigned int i = 1;	 // cpp11 does not support declaring the iterator in the
						 // lambda function
	std::generate(modes.begin(), modes.end(), [&]() mutable {
		i++;
		return f0 * i;
	});
	return modes;
}