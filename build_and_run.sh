#!/bin/bash
cmake -S . -B build
cmake -B build -DUSE_CLOUD=OFF
cmake --build build

./build/main
touch ./viewer/viewer.html

echo "finished"