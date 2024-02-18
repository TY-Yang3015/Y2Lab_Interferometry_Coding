cmake_minimum_required(VERSION 3.10)
project(InterferogramSimulation)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



add_executable(InterferogramSimulation main.cpp)


target_compile_options(InterferogramSimulation PRIVATE -O3)

# change this to your path when using!
include_directories(path/to/eigen-3.4.0)


