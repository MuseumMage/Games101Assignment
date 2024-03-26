cmake_minimum_required(VERSION 3.10)
project(Assignment4_BezierCurve)

set(CMAKE_CXX_STANDARD 14)

# OpenCV_DIR目录下要有OpenCVConfig.cmake，cmake才能找到opencv
set(OpenCV_DIR ${CMAKE_SOURCE_DIR}/lib/opencv)
#set(OpenCV_DIR ${CMAKE_SOURCE_DIR}/../lib/opencv)
# find_package 会从 (PACKAGENAME)_DIR去找
find_package(OpenCV REQUIRED)
message(STATUS "opencv include path is: " ${OpenCV_INCLUDE_DIRS})
message(STATUS "opencv lib path is: " ${OpenCV_LIBS})
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(Assignment4_BezierCurve main.cpp)

target_link_libraries(Assignment4_BezierCurve ${OpenCV_LIBRARIES})
