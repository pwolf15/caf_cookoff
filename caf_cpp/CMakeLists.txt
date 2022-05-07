cmake_minimum_required(VERSION 3.0)

project(caf_cpp)

add_executable(caf_cpp caf.cpp)
target_link_libraries(caf_cpp fftw3)