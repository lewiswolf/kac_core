cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 11)

# standalone
project(geometry_test)
add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/test.cpp)
