cmake_minimum_required(VERSION 3.10)

project(BendyPrint)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
find_package(LIBIGL REQUIRED QUIET)

set (CMAKE_C_STANDARD 99)

FILE(GLOB SRCFILES src/*.cpp)
FILE(GLOB HDRFILES src/*.h)

add_executable(${PROJECT_NAME} ${SRCFILES} ${HDRFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
