cmake_minimum_required(VERSION 3.16)
project(backmap_lambda_unit_tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
)
# Keep gtest's own runtime/CRT settings out of our way on all platforms.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(test_backmap_lambda test_backmap_lambda.cpp)
target_include_directories(test_backmap_lambda PRIVATE ${CMAKE_SOURCE_DIR}/../../src)
target_link_libraries(test_backmap_lambda PRIVATE gtest_main)

include(GoogleTest)
gtest_discover_tests(test_backmap_lambda)
