cmake_minimum_required(VERSION 3.12.0)
project(nbody LANGUAGES CXX C)

option(BUILD_SIMULATION "Build the simulation" ON)
option(BUILD_VISUALISATION "Build the simulation visualisation" OFF)
option(USE_OPENMP "Use OpenMP in the NBody solvers" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
message(STATUS)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_definitions(-DBUILD_TYPE_DEBUG)
endif()

if(BUILD_SIMULATION)

    # Find OpenMP
    if(USE_OPENMP)
        find_package(OpenMP)
        if(OpenMP_CXX_FOUND)
            message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
            message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")
        else()
            message(WARNING "OpenMP not found - building without parallel support")
            set(USE_OPENMP OFF)
        endif()
    endif()

    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

    add_executable(generate src/generate-input.cpp)
    add_executable(compareOutputs src/compare-outputs.cpp)
    add_executable(simulation src/simulation.cpp)

    # Link OpenMP if available
    if(USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(simulation OpenMP::OpenMP_CXX)
        target_compile_definitions(simulation PRIVATE USE_OPENMP=1)
        message(STATUS "Simulation built with OpenMP support")
    else()
        target_compile_definitions(simulation PRIVATE USE_OPENMP=0)
        message(STATUS "Simulation built WITHOUT OpenMP support")
    endif()

endif(BUILD_SIMULATION)

if(BUILD_VISUALISATION)
    add_subdirectory(external)

    find_package(OpenGL REQUIRED)

    file(GLOB SOURCES_VISUALISATION
        ${CMAKE_SOURCE_DIR}/src/vis/*.cpp
        ${CMAKE_SOURCE_DIR}/src/vis/shaders/*.cpp
        ${CMAKE_SOURCE_DIR}/src/vis/renderBufferObjects/*.cpp
        ${CMAKE_SOURCE_DIR}/src/vis/camera/*.cpp
    )

    add_executable(nbody_visualisation ${SOURCES_VISUALISATION})

    target_link_libraries(nbody_visualisation 
        glad
        glfw
        glm::glm
        OpenGL::GL
        imgui
        ${CMAKE_DL_LIBS}
    )

    add_custom_target(copy-shaders ALL 
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/src/vis/shaders/glsl" "${CMAKE_CURRENT_BINARY_DIR}/shaders" 
        COMMENT "Copying shaders"
    )
    add_dependencies(nbody_visualisation copy-shaders)

endif(BUILD_VISUALISATION)
