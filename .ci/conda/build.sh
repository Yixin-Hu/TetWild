#!/usr/bin/env bash
mkdir build
cd build

cmake -G "Ninja" \
	-DCMAKE_BUILD_TYPE="Release" \
	-DCMAKE_INSTALL_PREFIX:FILEPATH=$HOME/.local/ \
	-DTETWILD_WITH_ISPC=OFF \
	..

ninja install
