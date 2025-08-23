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
	std::string name;
	std::chrono::time_point<std::chrono::steady_clock> start_tp;
	// constructors
	Timer() { start_tp = std::chrono::steady_clock::now(); }
	explicit Timer(std::string s = ""): name(std::move(s)) {
		start_tp = std::chrono::steady_clock::now();
	}
	// destructors
	~Timer() {
		auto end_tp = std::chrono::steady_clock::now();
		auto duration = duration_cast<std::chrono::nanoseconds>(end_tp - start_tp).count();
		std::cout << (name.empty() ? "" : name + ": ") << duration / 1000 << "us\n";
	}
};

void booleanTest(const std::string& test_name, const bool& b) {
	/*
	Prints a warning to the console if the test fails.
	*/

	if (!b) {
		std::cerr << test_name << std::endl;
		throw;
	}
}

void batchBooleanTest(
	const std::string& test_name,
	const unsigned long& _N,
	const std::function<bool(const unsigned long&)>& lambda
) {
	/*
	Runs booleanTest N times on the input function.
	*/

	for (unsigned long _n = 0; _n < _N; _n++) {
		if (!lambda(_n)) {
			booleanTest(test_name, false);
			return;
		}
	}
}
