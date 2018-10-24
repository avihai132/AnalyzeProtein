cmake_minimum_required(VERSION 3.8)
project(AnalyzeProtein)

set(CMAKE_C_STANDARD 99)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wextra -Wall -Wvla -lm -std=c99")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

set(SOURCE_FILES AnalyzeProtein.c)
add_executable(AnalyzeProtein ${SOURCE_FILES})