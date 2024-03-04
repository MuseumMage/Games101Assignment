cmake_minimum_required (VERSION 2.8.11)
project (Transformation)

add_executable (Transformation main.cpp)

# Debug Message
# message("The project source directory is: ${CMAKE_SOURCE_DIR}/../lib")

target_include_directories(Transformation PUBLIC ${CMAKE_SOURCE_DIR}/../lib)