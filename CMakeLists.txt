cmake_minimum_required(VERSION 3.10)
project(OpenGL_Minimal)

set(CMAKE_CXX_STANDARD 17)

# Add glad
add_library(glad glad.c)
target_include_directories(glad PUBLIC include)

find_package(OpenGL REQUIRED)

file(GLOB SOURCES
    ${CMAKE_SOURCE_DIR}/*.cpp
    ${CMAKE_SOURCE_DIR}/shaders/*.cpp
    ${CMAKE_SOURCE_DIR}/particle/*.cpp
)

add_executable(nbody ${SOURCES})

target_link_libraries(nbody
    glad
    glfw
    OpenGL::GL
    ${CMAKE_DL_LIBS}
)
