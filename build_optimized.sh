#!/bin/bash

# Script to build FiniteFlow with performance optimizations

echo "Building FiniteFlow with performance optimizations..."

# Create build directory
mkdir -p build
cd build

# Configure with optimizations
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DFFLOW_THREAD_POOL=ON \
    -DBUILD_SHARED_LIBS=ON

# Build with all available cores
make -j$(nproc)

echo "Build completed! The optimized library is in the build directory."
echo "To install the library system-wide, run: sudo make install"