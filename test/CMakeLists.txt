cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 11)

project(geometry_test)

# include the test src dir
include_directories(
	"${CMAKE_CURRENT_SOURCE_DIR}/../src"
)
# include geometry.hpp
file(GLOB_RECURSE GEOMETRY
	"${CMAKE_CURRENT_SOURCE_DIR}/../src/geometry.hpp"
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/test.cpp)
