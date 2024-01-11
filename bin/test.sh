#!/bin/bash

# build dir
if [ ! -d ./build ]; then
  mkdir -p ./build
fi

# build
cmake -S . -B build
cmake --build build --config Release -j

# run test
ctest --test-dir build --config Release -j --output-on-failure