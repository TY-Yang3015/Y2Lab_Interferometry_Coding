cmake_minimum_required(VERSION 3.24)
project(CrossingPoints)

set(CMAKE_CXX_STANDARD 17)


#change this to your path when using
include_directories(path/to/eigen-3.4.0)



add_executable(CrossingPoints main.cpp)

find_package(iir)
target_link_libraries(CrossingPoints iir::iir)


target_compile_options(CrossingPoints PRIVATE -O3)