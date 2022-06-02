cmake -S . -B build
cmake --build build
cmake -E env CTEST_OUTPUT_ON_FAILURE=1 cmake --build build --target test