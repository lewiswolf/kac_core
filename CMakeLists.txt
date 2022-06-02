cmake_minimum_required(VERSION 3.22)
project(
	kac_core
	VERSION 0.0.2
	DESCRIPTION ""
	LANGUAGES CXX
)

# create library
add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

# link /src
target_include_directories(${PROJECT_NAME} INTERFACE src)
target_sources(
	${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/kac_core.hpp
)

# include dependencies
include(FetchContent)
FetchContent_Declare(
  boost_math
  GIT_REPOSITORY https://github.com/boostorg/math.git
)
FetchContent_MakeAvailable(boost_math)
target_link_libraries(${PROJECT_NAME} INTERFACE boost_math)

# if this is project root, run tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	include(CTest)
	set(CMAKE_CXX_STANDARD 11)
	add_subdirectory(test)
endif()