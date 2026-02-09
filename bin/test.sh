#!/bin/bash
# Run build scripts.
set -e

# Parse arguments
BUILD_TYPE="Debug"
CTEST_ARGS=""
while getopts "DRV" opt; do
	case "$opt" in
		D)
			BUILD_TYPE="Debug"
			;;
		R)
			BUILD_TYPE="Release"
			;;
		V)
			CTEST_ARGS="--verbose --output-on-failure"
			;;
		*)
			echo "Usage: $0 [-D] [-R] [-V]"
			echo "-D causes ctest to build the library in Debug (default)."
			echo "-R causes ctest to build the library in Release."
			echo "-V causes ctest to print it's the entire output including calls to std::cout."
			exit 1
			;;
	esac
done

# build dir
if [ ! -d ./build ]; then
	mkdir -p ./build
fi

# build
cmake -S . -B build -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
cmake --build build -j

# run test
ctest --test-dir build --build-config "$BUILD_TYPE" -j $CTEST_ARGS
echo
./build/test/profiler