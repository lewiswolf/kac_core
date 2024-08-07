#!/bin/bash

# build dir
if [ ! -d ./build ]; then
	mkdir -p ./build
fi

# build
cmake -S . -B build
cmake --build build --config Debug -j

# run test
./build/test/profiler
echo 
ctest --test-dir build --build-config Debug -j --output-on-failure