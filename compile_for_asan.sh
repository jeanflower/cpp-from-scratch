#!/bin/bash

# Navigate to the build directory
cd build
echo "Clean the build directory (remove all files in the directory)"
rm -rf *
cd ..

echo "Build with clang"
#clang++ -std=c++17 -fsanitize=address -g -o build/main main.cpp
clang++ -std=c++17 -g -o build/main \
    -arch x86_64 \
    -I/usr/local/Cellar/opencascade/7.8.1_1/include/opencascade \
    -L/usr/local/Cellar/opencascade/7.8.1_1/lib \
    -lTKernel -lTKMath -lTKG2d -lTKG3d -lTKGeomBase -lTKGeomAlgo \
    main.cpp \
    boost_explore/boost_usage.cpp \
    std_explore/std_usage.cpp \
    geom_explore/geomAPI_usage.cpp

# Run the program and analyze leaks
leaks -atExit -- ./build/main

#Look for output like this if you're looking for a leak
#
#Process 58508: 1 leak for 16 total leaked bytes.
#
#STACK OF 1 INSTANCE OF 'ROOT LEAK: <malloc in boost_data_types::shared_ptr_example()>':
#4   dyld                                  0x19310b154 start + 2476
#3   main                                  0x100a046e8 main + 32  main.cpp:5
#2   main                                  0x100a04710 boost_data_types::shared_ptr_example() + 24  boost_usage.cpp:10
#1   libc++abi.dylib                       0x19344ebd4 operator new(unsigned long) + 32
#0   libsystem_malloc.dylib                0x1932cfa68 _malloc_zone_malloc_instrumented_or_legacy + 148 
#====
#    1 (16 bytes) ROOT LEAK: <malloc in boost_data_types::shared_ptr_example() 0x12b7040f0> [16]
