#!/usr/bin/env bash
mkdir build
cd build

cmake -G "Ninja" \
    -DCMAKE_BUILD_TYPE="Release" \
	-DCMAKE_INSTALL_PREFIX:FILEPATH=$HOME/.local/ \
	-DBUILD_TESTING=OFF \
    -DUSE_ISPC=OFF \
	..

ninja install
