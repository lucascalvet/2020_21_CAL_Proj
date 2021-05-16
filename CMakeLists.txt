cmake_minimum_required(VERSION 3.19)
project(2020_21_CAL_Proj)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(2020_21_CAL_Proj main.cpp Graph.h MutablePriorityQueue.h)
