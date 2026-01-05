cmake_minimum_required(VERSION 3.20.1)
project(test-ampere)
set(SOURCES test_ampere.cpp
    ${CMAKE_SOURCE_DIR}/src/vecfield.hpp
    ${CMAKE_SOURCE_DIR}/src/field.hpp
    ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
    ${CMAKE_SOURCE_DIR}/src/population.hpp
    ${CMAKE_SOURCE_DIR}/src/ampere.hpp
)
add_executable(${PROJECT_NAME} ${SOURCES})
add_test(NAME test-ampere COMMAND test-ampere)
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
message(${PROJECT_NAME} " target: ${SOURCES}")
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
