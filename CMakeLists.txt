################################################################################
# TetWild CMake
################################################################################
cmake_minimum_required(VERSION 3.3)
project(TetWild)

################################################################################
# Options
################################################################################
# tetwild
option(TETWILD_QUIET     "Mute log and unnecassary checks" OFF)
option(TETWILD_WITH_ISPC "Use ISPC"                        OFF)
# libigl library
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR      "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           ON)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK             "Use Cork"           OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          OFF)
option(LIBIGL_WITH_PNG              "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON           "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       OFF)
option(LIBIGL_WITH_XML              "Use XML"            OFF)
# gui
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  OFF)
#geogram
option(GEOGRAM_WITH_TRIANGLE        "Use Triangle"       OFF)

################################################################################
# Settings
################################################################################
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(TETWILD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/extern")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# 3rd party libraries
################################################################################
# cgal
find_package(CGAL REQUIRED)

# libigl
find_package(LIBIGL REQUIRED)

# geogram
include(geogram)

# pymesh loaders
add_subdirectory(${TETWILD_EXTERNAL}/pymesh)

# CL11
add_subdirectory(${TETWILD_EXTERNAL}/CLI)

# cout
if(TETWILD_QUIET)
	add_definitions(-DMUTE_COUT)
endif()

################################################################################
# TetWild
################################################################################
# Build static library for executable
add_library(libTetWild STATIC
		src/tetwild/Preprocess.h
		src/tetwild/Preprocess.cpp
		src/tetwild/DelaunayTetrahedralization.h
		src/tetwild/DelaunayTetrahedralization.cpp
		src/tetwild/MeshConformer.h
		src/tetwild/MeshConformer.cpp
		src/tetwild/BSPSubdivision.h
		src/tetwild/BSPSubdivision.cpp
		src/tetwild/SimpleTetrahedralization.h
		src/tetwild/SimpleTetrahedralization.cpp
		src/tetwild/MeshRefinement.h
		src/tetwild/MeshRefinement.cpp
		src/tetwild/LocalOperations.h
		src/tetwild/LocalOperations.cpp
		src/tetwild/EdgeSplitter.h
		src/tetwild/EdgeSplitter.cpp
		src/tetwild/EdgeCollapser.h
		src/tetwild/EdgeCollapser.cpp
		src/tetwild/EdgeRemover.h
		src/tetwild/EdgeRemover.cpp
		src/tetwild/VertexSmoother.h
		src/tetwild/VertexSmoother.cpp
		src/tetwild/InoutFiltering.h
		src/tetwild/InoutFiltering.cpp
		src/tetwild/heads.h
		src/tetwild/heads.cpp
		src/tetwild/tetwild.cpp
		src/tetwild/TetmeshElements.h
		include/tetwild/tetwild.h
)
target_include_directories(libTetWild
	PRIVATE
		src
	PUBLIC
		include
)
target_link_libraries(libTetWild
	PUBLIC
		geogram
		igl::core
		igl::cgal
		pymesh::pymesh
)
set_target_properties(libTetWild PROPERTIES OUTPUT_NAME "tetwild")

# ispc
if(TETWILD_WITH_ISPC)
	message(STATUS "Compiling energy with ISPC")
	add_subdirectory(src/ispc)
	ispc_add_energy(libTetWild)
endif()

# Building executable
add_executable(TetWild src/main.cpp)
target_link_libraries(TetWild
	libTetWild
	CLI11::CLI11
)
target_include_directories(TetWild PRIVATE src)

# Install
install(TARGETS TetWild RUNTIME DESTINATION bin)
