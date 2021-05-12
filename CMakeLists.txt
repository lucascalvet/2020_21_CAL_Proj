cmake_minimum_required(VERSION 3.19)
project(2020_21_CAL_Proj)

set(CMAKE_C_STANDARD 11)

# Add the source files
file(GLOB SRC_FILES CONFIGURE_DEPENDS "*.cpp" "*.h")
add_executable(2020_21_CAL_Proj ${SRC_FILES})
