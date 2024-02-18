cmake_minimum_required(VERSION 3.10)
project(InterferogramSimulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Add the executable
add_executable(InterferogramSimulation main.cpp)

# Example of adding optimization flags and other compiler options
target_compile_options(InterferogramSimulation PRIVATE -O3)

include_directories(/Users/yangtianyue/Desktop/Programming/CPP_Libs/eigen-3.4.0)


