cmake_minimum_required (VERSION 2.8.11)
project (Assignment0_Transformation)

add_executable (Assignment0_Transformation main.cpp)

# Debug Message
# message("The project source directory is: ${CMAKE_SOURCE_DIR}/../lib")

target_include_directories(Assignment0_Transformation PUBLIC ${CMAKE_SOURCE_DIR}/../lib)