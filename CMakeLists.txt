cmake_minimum_required(VERSION 3.1)

project(single_e)

# boost
find_package(Boost REQUIRED)

# python
find_package(PythonLibs REQUIRED)

# Pybind11
add_subdirectory(pybind11)

pybind11_add_module(single_e single_e.cpp)

target_link_libraries(single_e PRIVATE "${Boost_LIBRARIES}")
target_include_directories(single_e PRIVATE ${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})

add_custom_command(
    TARGET single_e POST_BUILD VERBATIM
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:single_e>" "${CMAKE_SOURCE_DIR}/isolver.so"
)
