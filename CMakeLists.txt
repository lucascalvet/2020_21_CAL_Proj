cmake_minimum_required(VERSION 3.19)
project(2020_21_CAL_Proj)

set(CMAKE_C_STANDARD 11)

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Build GraphViewerCpp
set(GRAPHVIEWERCPP_PATH lib/GraphViewerCpp)
add_subdirectory(${GRAPHVIEWERCPP_PATH})
include_directories(${GRAPHVIEWERCPP_PATH}/include)

include_directories(.)

# Add the source files
file(GLOB SRC_FILES CONFIGURE_DEPENDS "*.cpp" "*.h")
add_executable(2020_21_CAL_Proj main.cpp Graph.h Utils.cpp Utils.h Menu.h Menu.cpp)
target_link_libraries(2020_21_CAL_Proj graphviewer)
