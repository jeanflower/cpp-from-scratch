cmake_minimum_required(VERSION 3.12)

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

# Set the policy to handle package root variables in uppercase
cmake_policy(SET CMP0144 NEW)
cmake_policy(SET CMP0167 NEW)

# Find Boost package
set(BOOST_ROOT "/usr/local/Cellar/boost/1.87.0/")
find_package(Boost 1.70 REQUIRED COMPONENTS filesystem system)

# Manually include Catch2 headers
include_directories("/usr/local/Cellar/catch2/3.8.0/include")

# Add executable target (your main application)
add_executable(main main.cpp)

# Add include directories
include_directories(${Boost_INCLUDE_DIRS})

# Link Boost libraries to your project
target_link_libraries(main ${Boost_LIBRARIES})

# Add a test executable (your test file)
add_executable(test_main test_main.cpp)

# Link Catch2 with the test executable
target_link_libraries(test_main Catch2::Catch2WithMain)

# Add a test for CTest
add_test(NAME MyTest COMMAND test_main)
