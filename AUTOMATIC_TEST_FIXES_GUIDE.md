# Automatic Test Fixes Guide

This document shows the exact changes needed for each PR to add CTest integration.

## Changes Required for ALL PRs

### 1. Update Root CMakeLists.txt

The root `CMakeLists.txt` currently has:
```cmake
enable_testing()
add_subdirectory(tests/boris)
```

It should be updated to include ALL test subdirectories that exist in each PR:

```cmake
enable_testing()
add_subdirectory(tests/boris)
add_subdirectory(tests/ampere)           # If exists in PR
add_subdirectory(tests/faraday)          # If exists in PR
add_subdirectory(tests/population)       # If exists in PR
add_subdirectory(tests/ampere_faraday)   # If exists in PR
```

### 2. Add `add_test()` Commands to Each Test CMakeLists.txt

## PR #3 (Tomas Formanek) - AUTOMATIC

### tests/boris/CMakeLists.txt

**Current:**
```cmake
cmake_minimum_required(VERSION 3.20.1)
project(test-boris)
set(SOURCES test_boris.cpp
    ${CMAKE_SOURCE_DIR}/src/pusher.hpp
    ${CMAKE_SOURCE_DIR}/src/vecfield.hpp
    ${CMAKE_SOURCE_DIR}/src/field.hpp
    ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
)
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
message(${PROJECT_NAME} " target: ${SOURCES}")
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
```

**Should Be:**
```cmake
cmake_minimum_required(VERSION 3.20.1)
project(test-boris)
set(SOURCES test_boris.cpp
    ${CMAKE_SOURCE_DIR}/src/pusher.hpp
    ${CMAKE_SOURCE_DIR}/src/vecfield.hpp
    ${CMAKE_SOURCE_DIR}/src/field.hpp
    ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
)
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)

# Add CTest integration
add_test(NAME boris_pusher_test COMMAND ${PROJECT_NAME})
```

### tests/ampere_faraday/CMakeLists.txt

**Current:**
```cmake
cmake_minimum_required(VERSION 3.20.1)
project(test-ampere-faraday)
set(SOURCES test_ampere_faraday.cpp
    ${CMAKE_SOURCE_DIR}/src/ampere.hpp
    ${CMAKE_SOURCE_DIR}/src/faraday.hpp
    ${CMAKE_SOURCE_DIR}/src/vecfield.hpp
    ${CMAKE_SOURCE_DIR}/src/field.hpp
    ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
    ${CMAKE_SOURCE_DIR}/src/boundary_condition.hpp
)
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
message(${PROJECT_NAME} " target: ${SOURCES}")
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
```

**Should Be:**
```cmake
cmake_minimum_required(VERSION 3.20.1)
project(test-ampere-faraday)
set(SOURCES test_ampere_faraday.cpp
    ${CMAKE_SOURCE_DIR}/src/ampere.hpp
    ${CMAKE_SOURCE_DIR}/src/faraday.hpp
    ${CMAKE_SOURCE_DIR}/src/vecfield.hpp
    ${CMAKE_SOURCE_DIR}/src/field.hpp
    ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
    ${CMAKE_SOURCE_DIR}/src/boundary_condition.hpp
)
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)

# Add CTest integration
add_test(NAME ampere_faraday_test COMMAND ${PROJECT_NAME})
```

### Root CMakeLists.txt

Add after `enable_testing()`:
```cmake
add_subdirectory(tests/ampere_faraday)
```

---

## PR #15 (Gaia Bilosi) - AUTOMATIC

This PR has FOUR test directories, each needs `add_test()`:

### tests/boris/CMakeLists.txt

Add at end:
```cmake
add_test(NAME boris_pusher_test COMMAND ${PROJECT_NAME})
```

### tests/ampere/CMakeLists.txt

Add at end:
```cmake
add_test(NAME ampere_law_test COMMAND ${PROJECT_NAME})
```

### tests/faraday/CMakeLists.txt

Add at end:
```cmake
add_test(NAME faraday_law_test COMMAND ${PROJECT_NAME})
```

### tests/population/CMakeLists.txt

Add at end:
```cmake
add_test(NAME population_moments_test COMMAND ${PROJECT_NAME})
```

### Root CMakeLists.txt

Update to:
```cmake
enable_testing()
add_subdirectory(tests/boris)
add_subdirectory(tests/ampere)
add_subdirectory(tests/faraday)
add_subdirectory(tests/population)
```

---

## PR #4, #5, #6, #8, #10, #11, #12, #14 - AUTOMATIC

Each of these PRs needs similar fixes. The pattern is:

1. Find all test subdirectories in the PR
2. Add `add_test(NAME <descriptive-name> COMMAND ${PROJECT_NAME})` to each test's CMakeLists.txt
3. Add all test subdirectories to root CMakeLists.txt

---

## How to Verify the Fix Works

After adding `add_test()` commands:

```bash
mkdir build && cd build
cmake ..
make
ctest --verbose
```

Expected output:
```
Test project /path/to/build
    Start 1: boris_pusher_test
1/4 Test #1: boris_pusher_test ................   Passed    0.05 sec
    Start 2: ampere_law_test
2/4 Test #2: ampere_law_test ..................   Passed    0.03 sec
    Start 3: faraday_law_test
3/4 Test #3: faraday_law_test .................   Passed    0.03 sec
    Start 4: population_moments_test
4/4 Test #4: population_moments_test ..........   Passed    0.04 sec

100% tests passed, 0 tests failed out of 4
```

---

## GitHub Actions Integration

To make tests run automatically in CI, ensure `.github/workflows/*.yml` includes:

```yaml
- name: Run tests
  run: |
    cd build
    ctest --output-on-failure
```

This will cause GitHub Actions to fail if any test fails, ensuring code quality.

---

## Summary

**For EACH PR, create a new PR titled:**
- "PR #3 Tomas Formanek AUTOMATIC"
- "PR #4 Boris pusher AUTOMATIC"
- "PR #15 HPC PROJECT - Bilosi Gaia AUTOMATIC"
- etc.

**Each AUTOMATIC PR should:**
1. Add `add_test()` command to every test CMakeLists.txt
2. Update root CMakeLists.txt to add_subdirectory for all tests
3. Include a brief description explaining the CTest integration

**Result:** Tests will run automatically in GitHub Actions CI pipeline.
