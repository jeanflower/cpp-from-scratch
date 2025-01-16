#!/bin/bash

# Make sure you're in the correct build directory
cd build

# If necessary, rebuild the project (optional)
# cmake ..
# cmake --build .

# Start the debugger (lldb for macOS or gdb for Linux)
# For macOS:
lldb ./main

# For Linux (use gdb instead):
# gdb ./main
