# tests need to be added as executables first
add_executable(profiler src/profiler.cpp)
add_executable(test_fdtd src/test_fdtd.cpp)
add_executable(test_geometry src/test_geometry.cpp)
add_executable(test_modes src/test_modes.cpp)

# should be linked to the main library
target_link_libraries(profiler PRIVATE kac_core)
target_link_libraries(test_fdtd PRIVATE kac_core)
target_link_libraries(test_geometry PRIVATE kac_core)
target_link_libraries(test_modes PRIVATE kac_core)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME test_fdtd COMMAND test_fdtd)
add_test(NAME test_geometry COMMAND test_geometry)
add_test(NAME test_modes COMMAND test_modes)