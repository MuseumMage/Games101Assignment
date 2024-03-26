cmake_minimum_required(VERSION 3.10)
project(Assignment1_Rasterizer)
set(CMAKE_CXX_STANDARD 17)

# OpenCV_DIR目录下要有OpenCVConfig.cmake，cmake才能找到opencv
set(OpenCV_DIR ${CMAKE_SOURCE_DIR}/lib/opencv)
message("zhangbo" ${OpenCV_DIR})
#set(OpenCV_DIR ${CMAKE_SOURCE_DIR}/../lib/opencv)
# find_package 会从 (PACKAGENAME)_DIR去找
find_package(OpenCV REQUIRED)

message(STATUS "opencv include path is: " ${OpenCV_INCLUDE_DIRS})
message(STATUS "opencv lib path is: " ${OpenCV_LIBS})
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(Assignment1_Rasterizer main.cpp rasterizer.hpp rasterizer.cpp Triangle.hpp Triangle.cpp)
target_link_libraries(Assignment1_Rasterizer ${OpenCV_LIBS})
target_include_directories(Assignment1_Rasterizer PUBLIC ${CMAKE_SOURCE_DIR}/../lib)
