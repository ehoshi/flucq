#!/bin/sh

mysrcdir=`pwd`
mybuilddir=`pwd`/build
FC=gfortran
export FC
mkdir -p "$mybuilddir"
cd "$mybuilddir"
cmake -D CMAKE_BUILD_TYPE=Release "$mysrcdir"
make flucq
make props
cp flucq flucq_Linux
cp props props_Linux
