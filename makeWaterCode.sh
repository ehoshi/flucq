#!/bin/bash

mysrcdir=`pwd`
MD_ARCH=`uname -s`
export MD_ARCH
FC=gfortran
export FC
make flucq
make props
