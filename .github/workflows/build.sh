#!/bin/bash

set -e

mkdir build
cd build
cmake ..
make -j2
