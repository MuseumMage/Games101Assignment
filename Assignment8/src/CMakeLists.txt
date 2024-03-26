cmake_minimum_required(VERSION 3.10)

# Application source
set(APPLICATION_SOURCE
    rope.cpp
    application.cpp
    main.cpp
    getopt.c
    getopt_init.c
    ./../CGL/deps/glew/src/glew.c
    ./../CGL/src/osdfont.c
)

#-------------------------------------------------------------------------------
# Set include directories
#-------------------------------------------------------------------------------
include_directories(
  ${GLEW_INCLUDE_DIRS}
  ${GLFW_INCLUDE_DIRS}
  ${CGL_INCLUDE_DIRS}
  ${FREETYPE_INCLUDE_DIRS}
)

#-------------------------------------------------------------------------------
# Set link directories
#-------------------------------------------------------------------------------
link_directories(
  ${GLEW_LIBRARY_DIRS}
  ${GLFW_LIBRARY_DIRS}
  ${CGL_LIBRARY_DIRS}
  ${FREETYPE_LIBRARY_DIRS}
)

# test
message("zhangbo GLEW_INCLUDE_DIRS: ${GLEW_INCLUDE_DIRS}")
message("zhangbo GLFW_INCLUDE_DIRS: ${GLFW_INCLUDE_DIRS}")
message("zhangbo CGL_INCLUDE_DIRS: ${CGL_INCLUDE_DIRS}")
message("zhangbo GLEW_LIBRARY_DIRS: ${GLEW_LIBRARY_DIRS}")
message("zhangbo GLFW_LIBRARY_DIRS: ${GLFW_LIBRARY_DIRS}")
message("zhangbo CGL_LIBRARY_DIRS: ${CGL_LIBRARY_DIRS}")

#-------------------------------------------------------------------------------
# Add executable
#-------------------------------------------------------------------------------
add_executable(Assignment8_RopeSim ${APPLICATION_SOURCE})

target_link_libraries( Assignment8_RopeSim
    CGL ${CGL_LIBRARIES}
    glew ${GLEW_LIBRARIES}
    glfw ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${FREETYPE_LIBRARIES}
    ${CMAKE_THREADS_INIT}
    freetype
)

#-------------------------------------------------------------------------------
# Platform-specific configurations for target
#-------------------------------------------------------------------------------
if(APPLE)
  set_property( TARGET ropesim APPEND_STRING PROPERTY COMPILE_FLAGS
                "-Wno-deprecated-declarations -Wno-c++11-extensions")
endif(APPLE)

# Put executable in build directory root
set(EXECUTABLE_OUTPUT_PATH ..)

# Install to project root
install(TARGETS Assignment8_RopeSim DESTINATION ${RopeSim_SOURCE_DIR})
