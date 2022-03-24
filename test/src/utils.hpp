#pragma once

// core
#include <chrono>
#include <iostream>
#include <string>

void booleanTest(const bool& b) {
	if (b) {
		std::cout << "✅\n";
	} else {
		std::cout << "❌\n";
		throw;
	}
}

class Timer {
   public:
	Timer() { start_tp = std::chrono::high_resolution_clock::now(); };
	Timer(std::string s) {
		name = s;
		start_tp = std::chrono::high_resolution_clock::now();
	}

	~Timer() { stop(); };

	void stop() {
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

   private:
	std::string name = "";
	std::chrono::time_point<std::chrono::high_resolution_clock> start_tp;
};