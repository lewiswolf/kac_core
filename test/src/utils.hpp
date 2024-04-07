/*
Utility functions for testing.
*/

#pragma once
// core
#include <chrono>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

struct Timer {
	/*
	Scoped timer, operating for as long as the scope is active. When
	the timer is destructed, it will print its lifetime to the console.
	*/

	// vars
	std::string name = "";
	std::chrono::time_point<std::chrono::high_resolution_clock> start_tp;
	// constructors
	Timer() { start_tp = std::chrono::high_resolution_clock::now(); }
	Timer(std::string s) {
		name = s;
		start_tp = std::chrono::high_resolution_clock::now();
	}
	// destructors
	~Timer() {
		auto end_tp = std::chrono::high_resolution_clock::now();
		auto start = std::chrono::time_point_cast<std::chrono::microseconds>(start_tp)
						 .time_since_epoch()
						 .count();
		auto end = std::chrono::time_point_cast<std::chrono::microseconds>(end_tp)
					   .time_since_epoch()
					   .count();
		std::cout << (name != "" ? name + ": " : "") << end - start << "us\n";
	}
};

void booleanTest(const std::string& test_name, const bool& b) {
	/*
	Prints a warning to the console if the test fails.
	*/

	if (!b) {
		std::cout << test_name;
		throw;
	}
}

void batchBooleanTest(
	const std::string& test_name,
	const unsigned long& _N,
	const std::function<bool(const unsigned long)>& lambda
) {
	/*
	Runs booleanTest N times on the input function.
	*/

	for (unsigned long n = 0; n < _N; n++) {
		if (!lambda(n)) {
			booleanTest(test_name, false);
			return;
		}
	}
	booleanTest(test_name, true);
}
