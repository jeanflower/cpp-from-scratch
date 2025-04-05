# From
# https://cloud.google.com/run/docs/quickstarts/build-and-deploy/deploy-c-plus-plus-service?hl=en
cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-std=c++23)

# creates profiling info which can be used by Instruments
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")

# Project name and version
project(MyProject VERSION 1.0)

# Set the policy to handle package root variables in uppercase
cmake_policy(SET CMP0144 NEW)
cmake_policy(SET CMP0167 NEW)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Add debugging and AddressSanitizer flags for debug builds
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-D_GLIBCXX_DEBUG -fsanitize=address -g)
  add_link_options(-fsanitize=address)
endif()

# Find Catch2 package
find_package(Catch2 3 REQUIRED)
# Find Boost package
set(BOOST_ROOT "/usr/local/Cellar/boost/1.87.0/")
find_package(Boost 1.70 REQUIRED COMPONENTS filesystem system)
find_package(OpenCASCADE REQUIRED)

# Add include directories
# Manually include Catch2 headers
include_directories("/usr/local/Cellar/catch2/3.8.0/include")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${OpenCASCADE_INCLUDE_DIR})

# Add the include directory for your headers
include_directories(boost_explore)
include_directories(std_explore)
include_directories(geom_explore)
include_directories(geomAPI_explore)
include_directories(nr_explore)
include_directories(rapidjson)
    
link_directories(${OpenCASCADE_LIBRARY_DIR})

# Add the source files
add_library(boost_explore_lib boost_explore/boost_usage.cpp)
add_library(std_explore_lib std_explore/std_usage.cpp)
add_library(geomAPI_lib geom_explore/geomAPI_usage.cpp)
add_library(geom_lib 
  geom_explore/geom.cpp
  geom_explore/geom_display.cpp
  geom_explore/geom_nurbs.cpp      
)
add_library(nr_lib nr_explore/gaussj.cpp)

# Add executable target (your main file)
add_executable(main main.cpp)

# Add a test executable (your test files)
add_executable(test_main test_main.cpp)
add_executable(test_geom test_geom.cpp)

# Link libraries
target_link_libraries(geom_lib 
    PRIVATE nr_lib
)

target_link_libraries(main 
    PRIVATE ${Boost_LIBRARIES} 
            TKernel TKMath TKGeomBase TKGeomAlgo
            boost_explore_lib 
            std_explore_lib
            geomAPI_lib
            geom_lib
            nr_lib
)

target_link_libraries(test_main 
  PRIVATE Catch2::Catch2WithMain
            ${Boost_LIBRARIES} 
            TKernel TKMath TKGeomBase TKGeomAlgo
            boost_explore_lib 
            std_explore_lib
            geomAPI_lib
            geom_lib
)

target_link_libraries(test_geom 
  PRIVATE Catch2::Catch2WithMain
          ${Boost_LIBRARIES} 
          geom_lib
          nr_lib
)

# Enable testing
enable_testing()
# Add a test for CTest
add_test(NAME MainTest COMMAND test_main)
add_test(NAME GeomTest COMMAND test_geom)
