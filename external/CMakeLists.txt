add_library(glad glad/glad.c)
target_include_directories(glad PUBLIC glad/include)

option(GLFW_LIBRARY TYPE STATIC)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_BUILD_DOCS OFF)
add_subdirectory(glfw)

include_directories(
    glfw/include
    glfw/deps
    )

add_subdirectory(glm)
