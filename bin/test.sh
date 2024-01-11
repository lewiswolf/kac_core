#!/bin/bash

# build dir
if [ ! -d ./build ]; then
  mkdir -p ./build
fi

# build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug -j

# run test
ctest --test-dir build --C Debug -j --output-on-failure