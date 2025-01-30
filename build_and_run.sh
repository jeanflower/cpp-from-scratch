#!/bin/bash
cd build
cmake -S .. -B build_local
cmake -B build_local -DUSE_CLOUD=OFF
cmake --build build_local

cd ..

./build/main
touch ./viewer/viewer.html

echo "finished"