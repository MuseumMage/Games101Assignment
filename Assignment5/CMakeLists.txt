cmake_minimum_required(VERSION 3.10)
project(Assignment5_RayTracing)

set(CMAKE_CXX_STANDARD 17)

add_executable(Assignment5_RayTracing main.cpp Object.hpp Vector.hpp Sphere.hpp global.hpp Triangle.hpp Scene.cpp Scene.hpp Light.hpp Renderer.cpp)
# target_compile_options(RayTracing PUBLIC -Wall -Wextra -pedantic -Wshadow -Wreturn-type -fsanitize=undefined)
#target_compile_options(RayTracing PUBLIC -Wall -pedantic -fsanitize=undefined)
# if(MSVC)
#    target_compile_options(RayTracing PUBLIC /W4 /WX)
# else()
#    target_compile_options(RayTracing PUBLIC -Wall -Wextra -pedantic -Wshadow -Wreturn-type -fsanitize=undefined)
# endif()

target_compile_features(Assignment5_RayTracing PUBLIC cxx_std_17)
target_link_libraries(Assignment5_RayTracing PUBLIC -fsanitize=undefined)
