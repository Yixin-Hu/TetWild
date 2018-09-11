mkdir build
cd build

cmake -G "%GENERATOR%" ^
	-DCMAKE_BUILD_TYPE="Release" ^
	-DCMAKE_PREFIX_PATH:FILEPATH="%PREFIX%" ^
	-DCMAKE_INSTALL_PREFIX:FILEPATH="%LIBRARY_PREFIX%" ^
	-DTETWILD_WITH_ISPC=OFF ^
	..

msbuild /m /p:Configuration=Release INSTALL.vcxproj
