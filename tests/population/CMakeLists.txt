add_executable(test_population test_population.cpp)

target_include_directories(test_population PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_population PRIVATE HighFive)

