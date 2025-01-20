cmake_minimum_required(VERSION 3.12)
# Set the policy to handle package root variables in uppercase
cmake_policy(SET CMP0144 NEW)
cmake_policy(SET CMP0167 NEW)

# Project name and version
project(MyProject VERSION 1.0)
# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Find Catch2 package
find_package(Catch2 3 REQUIRED)
# Find Boost package
set(BOOST_ROOT "/usr/local/Cellar/boost/1.87.0/")
find_package(Boost 1.70 REQUIRED COMPONENTS filesystem system)


# Manually include Catch2 headers
include_directories("/usr/local/Cellar/catch2/3.8.0/include")
# Add include directories
include_directories(${Boost_INCLUDE_DIRS})

# Add the include directory for your headers
include_directories(boost_explore)
include_directories(std_explore)

# Add the source files
add_library(boost_explore_lib boost_explore/boost_usage.cpp)
add_library(std_explore_lib std_explore/std_usage.cpp)

# Add executable target (your main application)
add_executable(main main.cpp)
# Add a test executable (your test file)
add_executable(test_main test_main.cpp)

# Link Boost libraries to your project
target_link_libraries(main ${Boost_LIBRARIES})
# Link the libraries to your main executable
target_link_libraries(main boost_explore_lib)
target_link_libraries(main std_explore_lib)
# Link Catch2 with the test executable
target_link_libraries(test_main Catch2::Catch2WithMain)

# Enable testing
enable_testing()
# Add a test for CTest
add_test(NAME MyTest COMMAND test_main)
