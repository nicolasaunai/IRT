# Testing Improvements for PR #10

This document describes the changes needed to enable CTest integration for the tests in PR #10.

## Summary

PR #10 ("HPCproject") added two new test subdirectories:
- `tests/population/`
- `tests/ampere_faraday/`

However, these tests were not properly registered with CTest, so they couldn't be discovered and run using the `ctest` command. Additionally, the `enable_testing()` call was missing from the root CMakeLists.txt.

## Changes Required

### 1. Enable CTest in Root CMakeLists.txt

Add `enable_testing()` call before the test subdirectories are added:

```cmake
target_link_libraries(hybirt PRIVATE HighFive)

enable_testing()

add_subdirectory(tests/boris)
```

**File:** `CMakeLists.txt` (line 76)

### 2. Add test registration for test-population

Add `add_test()` command after `add_executable()`:

```cmake
add_executable(${PROJECT_NAME} ${SOURCES})
add_test(NAME test-population COMMAND test-population)
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
```

**File:** `tests/population/CMakeLists.txt` (line 12)

### 3. Add test registration for test-ampere-faraday

Add `add_test()` command after `add_executable()`:

```cmake
add_executable(${PROJECT_NAME} ${SOURCES})
add_test(NAME test-ampere-faraday COMMAND test-ampere-faraday)
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
```

**File:** `tests/ampere_faraday/CMakeLists.txt` (line 11)

### 4. Add test registration for test-boris (consistency)

For consistency with the other tests, add `add_test()` command:

```cmake
add_executable(${PROJECT_NAME} ${SOURCES})
add_test(NAME test-boris COMMAND test-boris)
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
```

**File:** `tests/boris/CMakeLists.txt` (line 10)

### 5. Update .gitignore

Add the build directory to .gitignore:

```
build/
```

**File:** `.gitignore`

## Verification

After applying these changes:

1. Configure the project:
   ```bash
   mkdir build && cd build
   cmake ..
   ```

2. Build all targets:
   ```bash
   make -j4
   ```

3. List available tests:
   ```bash
   ctest -N
   ```
   
   Expected output:
   ```
   Test project /path/to/IRT/build
     Test #1: test-boris
     Test #2: test-population
     Test #3: test-ampere-faraday
   
   Total Tests: 3
   ```

4. Run all tests:
   ```bash
   ctest
   ```
   
   Expected output:
   ```
   Test project /path/to/IRT/build
       Start 1: test-boris
   1/3 Test #1: test-boris .......................   Passed    0.01 sec
       Start 2: test-population
   2/3 Test #2: test-population ..................   Passed    0.17 sec
       Start 3: test-ampere-faraday
   3/3 Test #3: test-ampere-faraday ..............   Passed    0.01 sec
   
   100% tests passed, 0 tests failed out of 3
   
   Total Test time (real) =   0.19 sec
   ```

## Integration with GitHub Actions

With these changes in place, GitHub Actions workflows can now run:

```yaml
- name: Configure
  run: cmake -B build
  
- name: Build
  run: cmake --build build
  
- name: Test
  run: ctest --test-dir build --output-on-failure
```

This will automatically discover and run all registered tests, providing feedback on test results in the PR.

## Complete Diff

The complete diff from PR #10 to include these changes:

```diff
diff --git a/.gitignore b/.gitignore
index ef707e4..154f4fa 100644
--- a/.gitignore
+++ b/.gitignore
@@ -1,2 +1,3 @@
 subprojects
 compile_commands.json
+build/
diff --git a/CMakeLists.txt b/CMakeLists.txt
index 69dbea0..47b3a0a 100644
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -73,6 +73,7 @@ add_executable(hybirt ${SOURCE_INC} ${SOURCE_CPP})
 
 target_link_libraries(hybirt PRIVATE HighFive)
 
+enable_testing()
 
 add_subdirectory(tests/boris)
 add_subdirectory(tests/population)
diff --git a/tests/ampere_faraday/CMakeLists.txt b/tests/ampere_faraday/CMakeLists.txt
index c61c22f..d193da1 100644
--- a/tests/ampere_faraday/CMakeLists.txt
+++ b/tests/ampere_faraday/CMakeLists.txt
@@ -9,6 +9,7 @@ set(SOURCES test_ampere_faraday.cpp
     ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
 )
 add_executable(${PROJECT_NAME} ${SOURCES})
+add_test(NAME test-ampere-faraday COMMAND test-ampere-faraday)
 target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
 message(${PROJECT_NAME} " target: ${SOURCES}")
 target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
diff --git a/tests/boris/CMakeLists.txt b/tests/boris/CMakeLists.txt
index 6e18c0f..edd3804 100644
--- a/tests/boris/CMakeLists.txt
+++ b/tests/boris/CMakeLists.txt
@@ -7,6 +7,7 @@ set(SOURCES test_boris.cpp
     ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
 )
 add_executable(${PROJECT_NAME} ${SOURCES})
+add_test(NAME test-boris COMMAND test-boris)
 target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
 message(${PROJECT_NAME} " target: ${SOURCES}")
 target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
diff --git a/tests/population/CMakeLists.txt b/tests/population/CMakeLists.txt
index ab8ed97..e596356 100644
--- a/tests/population/CMakeLists.txt
+++ b/tests/population/CMakeLists.txt
@@ -10,6 +10,7 @@ set(SOURCES test_population.cpp
     ${CMAKE_SOURCE_DIR}/src/population.hpp
 )
 add_executable(${PROJECT_NAME} ${SOURCES})
+add_test(NAME test-population COMMAND test-population)
 target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
 message(${PROJECT_NAME} " target: ${SOURCES}")
 target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
```

## Branch Information

These changes have been prepared on the `hpcproject-automatic` branch, which is based on PR #10.

To create a new PR with these changes:
1. The branch `hpcproject-automatic` contains all changes from PR #10 plus the testing improvements
2. Create a new PR from `hpcproject-automatic` to `master` with title: "HPCproject AUTOMATIC"
3. This PR will enable automated test execution via GitHub Actions
