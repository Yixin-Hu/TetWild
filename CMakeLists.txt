################################################################################
# TetWild CMake
################################################################################
cmake_minimum_required(VERSION 3.3)
project(TetWild)

################################################################################
# Options
################################################################################
# ispc
option(GTET_ISPC "Use ISPC" OFF)
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

################################################################################
# Settings
################################################################################
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

################################################################################
# 3rd party libraries
################################################################################
# cgal
find_package(CGAL REQUIRED)

# libigl
find_package(LIBIGL REQUIRED)

# geogram
include(geogram)

# ispc
SET(OBJ_FILES "")
if(GTET_ISPC)
	add_custom_command(OUTPUT energy_ispc.o
			COMMAND ispc ${CMAKE_SOURCE_DIR}/ispc/energy.ispc -h ${CMAKE_SOURCE_DIR}/ispc/energy.h -o energy_ispc.o DEPENDS ${CMAKE_SOURCE_DIR}/ispc/energy.ispc)
	add_definitions(-DGTET_ISPC)
    LIST(APPEND OBJ_FILES energy_ispc.o)
endif(GTET_ISPC)

################################################################################
# TetWild
################################################################################
include_directories(include)
# For PyMesh (change?)
include_directories(pymesh)

set(SOURCE_FILES
    src/Preprocess.h
	src/Preprocess.cpp
	src/DelaunayTetrahedralization.h
	src/DelaunayTetrahedralization.cpp
	src/MeshConformer.h
	src/MeshConformer.cpp
	src/BSPSubdivision.h
	src/BSPSubdivision.cpp
	src/SimpleTetrahedralization.h
	src/SimpleTetrahedralization.cpp
	src/MeshRefinement.h
	src/MeshRefinement.cpp
	src/LocalOperations.h
	src/LocalOperations.cpp
	src/EdgeSplitter.h
	src/EdgeSplitter.cpp
	src/EdgeCollapser.h
	src/EdgeCollapser.cpp
	src/EdgeRemover.h
	src/EdgeRemover.cpp
	src/VertexSmoother.h
	src/VertexSmoother.cpp
	src/InoutFiltering.h
	src/InoutFiltering.cpp
	src/heads.h
	src/heads.cpp
	src/tetwild.cpp
	src/CLI11.hpp
	src/TetmeshElements.h
	# PyMesh
	pymesh/MshLoader.h
	pymesh/MshLoader.cpp
	pymesh/MshSaver.h
	pymesh/MshSaver.cpp
	pymesh/Exception.h
    )

# Build static library for executable
add_library(libTetWild STATIC ${SOURCE_FILES} ${OBJ_FILES})
target_link_libraries(libTetWild geogram igl::core igl::cgal)
target_include_directories(libTetWild SYSTEM PUBLIC ${EIGEN_INCLUDE_DIRS})

# Building exectuable
set(MAIN_FILE src/main.cpp)
add_executable(TetWild ${MAIN_FILE})
target_link_libraries(TetWild libTetWild)

# Install
INSTALL(TARGETS TetWild RUNTIME DESTINATION bin)
