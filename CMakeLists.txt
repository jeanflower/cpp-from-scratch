cmake_minimum_required(VERSION 3.10)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Project name and version
project(MyProject VERSION 1.0)

# Enable testing
enable_testing()

# Find Catch2 package
find_package(Catch2 3 REQUIRED)

# Manually include Catch2 headers
include_directories("/usr/local/Cellar/catch2/3.8.0/include")

# Add executable target (your main application)
add_executable(main main.cpp)

# Add a test executable (your test file)
add_executable(test_main test_main.cpp)

# Link Catch2 with the test executable
target_link_libraries(test_main Catch2::Catch2WithMain)

# Add a test for CTest
add_test(NAME MyTest COMMAND test_main)
