cmake_minimum_required(VERSION 3.10)
project(Assignment2_Rasterizer)

set(CMAKE_CXX_STANDARD 17)

# OpenCV_DIR目录下要有OpenCVConfig.cmake，cmake才能找到opencv
set(OpenCV_DIR ${CMAKE_SOURCE_DIR}/lib/opencv)
#set(OpenCV_DIR ${CMAKE_SOURCE_DIR}/../lib/opencv)
# find_package 会从 (PACKAGENAME)_DIR去找
find_package(OpenCV REQUIRED)

message(STATUS "opencv include path is: " ${OpenCV_INCLUDE_DIRS})
message(STATUS "opencv lib path is: " ${OpenCV_LIBS})
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(Assignment2_Rasterizer main.cpp rasterizer.hpp rasterizer.cpp global.hpp Triangle.hpp Triangle.cpp)
target_link_libraries(Assignment2_Rasterizer ${OpenCV_LIBRARIES})
target_include_directories(Assignment2_Rasterizer PUBLIC ${CMAKE_SOURCE_DIR}/../lib)