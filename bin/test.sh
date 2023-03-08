#!/bin/bash

# build dir
if [ ! -d ./build ]; then
  mkdir -p ./build
fi

# build
cmake -S . -B build
cmake --build build -j

# run test
cmake -E env CTEST_OUTPUT_ON_FAILURE=1 cmake --build build --target test