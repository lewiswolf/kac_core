cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 11)

# standalone
project(geometry_test)
include_directories(
	"${CMAKE_CURRENT_SOURCE_DIR}/../src"
)

file(GLOB_RECURSE GEOMETRY
	"${CMAKE_CURRENT_SOURCE_DIR}/../src/geometry.hpp"
)
add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/test.cpp)
