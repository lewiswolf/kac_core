/*
Utility functions for testing.
*/

#pragma once
// core
#include <chrono>
#include <iostream>
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
	Timer() { start_tp = std::chrono::high_resolution_clock::now(); };
	Timer(std::string s) {
		name = s;
		start_tp = std::chrono::high_resolution_clock::now();
	}

	// denstructors
	~Timer() {
		auto end_tp = std::chrono::high_resolution_clock::now();
		auto start =
			std::chrono::time_point_cast<std::chrono::microseconds>(start_tp)
				.time_since_epoch()
				.count();
		auto end =
			std::chrono::time_point_cast<std::chrono::microseconds>(end_tp)
				.time_since_epoch()
				.count();
		std::cout << (name != "" ? name + ": " : "") << end - start << "us\n";
	};
};

void booleanTest(const bool& b) {
	/*
	Prints a tick or a cross to the console.
	*/

	b ? std::cout << "✅\n" : std::cout << "❌\n";
}

void batchBooleanTest(
	const int& I, const std::function<bool(unsigned int)>& lambda
) {
	/*
	Runs boolean test on a batch of functions.
	*/

	for (unsigned int i = 0; i < I; i++) {
		if (!lambda(i)) {
			booleanTest(false);
			return;
		}
	}
	booleanTest(true);
}