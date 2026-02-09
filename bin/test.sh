#!/bin/bash
# Run build scripts.
# -V causes ctest to print it's the entire output including calls to std::cout.

# build dir
if [ ! -d ./build ]; then
	mkdir -p ./build
fi

# build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j

# run test
./build/test/profiler
echo
if [ "$1" == "-V" ]; then
    ctest --test-dir build --build-config Debug -j --output-on-failure -V
else
    ctest --test-dir build --build-config Debug -j --output-on-failure
fi