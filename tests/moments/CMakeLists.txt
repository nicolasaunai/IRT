cmake_minimum_required(VERSION 3.20.1)
project(test-moments)
set(SOURCES test_moments.cpp
    ${CMAKE_SOURCE_DIR}/src/pusher.hpp
    ${CMAKE_SOURCE_DIR}/src/vecfield.hpp
    ${CMAKE_SOURCE_DIR}/src/field.hpp
    ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
    ${CMAKE_SOURCE_DIR}/src/population.hpp
    ${CMAKE_SOURCE_DIR}/src/moments.hpp
)
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
add_test(NAME test-moments COMMAND test-moments)
message(${PROJECT_NAME} " target: ${SOURCES}")
