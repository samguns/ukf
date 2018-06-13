cmake_minimum_required(VERSION 3.5)
project(ukf)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")

add_executable(ukf main.cpp ukf.cpp)