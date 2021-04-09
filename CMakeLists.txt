cmake_minimum_required(VERSION 3.10)
project(TCRMatch VERSION 1.0)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -g")
add_executable(tcrmatch src/tcrmatch.cpp)