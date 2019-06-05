cmake_minimum_required(VERSION 3.3)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

################################################################################


################################################################################
# Options
################################################################################
# tetwild
option(TETWILD_WITH_HUNTER "Use Hunter to download and configure Boost" OFF)
option(TETWILD_WITH_ISPC   "Use ISPC"                                   OFF)
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

if(TETWILD_WITH_HUNTER)
	# Needs to be set before the main project... argh =/
	include(HunterGate)
	HunterGate(
		URL "https://github.com/ruslo/hunter/archive/v0.23.25.tar.gz"
		SHA1 "cb75cce9a3a8d552e70e7118f3203eb4ac05c201"
	)
endif()

################################################################################
# Project name
################################################################################

project(TetWild)

################################################################################
# Settings
################################################################################
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(TETWILD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/extern")

# Color output
include(UseColors)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Extra warnings
include(Warnings)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# 3rd party libraries
################################################################################
include(TetWildDependencies)

################################################################################
# TetWild
################################################################################
# Build static library for executable
add_library(libTetWild 
		include/tetwild/Args.h
		include/tetwild/Exception.h
		include/tetwild/Logger.h
		include/tetwild/tetwild.h
		src/tetwild/BSPSubdivision.cpp
		src/tetwild/BSPSubdivision.h
		src/tetwild/CGALTypes.h
		src/tetwild/Common.cpp
		src/tetwild/Common.h
		src/tetwild/DelaunayTetrahedralization.cpp
		src/tetwild/DelaunayTetrahedralization.h
		src/tetwild/DistanceQuery.cpp
		src/tetwild/DistanceQuery.h
		src/tetwild/EdgeCollapser.cpp
		src/tetwild/EdgeCollapser.h
		src/tetwild/EdgeRemover.cpp
		src/tetwild/EdgeRemover.h
		src/tetwild/EdgeSplitter.cpp
		src/tetwild/EdgeSplitter.h
		src/tetwild/ForwardDecls.h
		src/tetwild/InoutFiltering.cpp
		src/tetwild/InoutFiltering.h
		src/tetwild/LocalOperations.cpp
		src/tetwild/LocalOperations.h
		src/tetwild/Logger.cpp
		src/tetwild/MeshConformer.cpp
		src/tetwild/MeshConformer.h
		src/tetwild/MeshRefinement.cpp
		src/tetwild/MeshRefinement.h
		src/tetwild/Preprocess.cpp
		src/tetwild/Preprocess.h
		src/tetwild/SimpleTetrahedralization.cpp
		src/tetwild/SimpleTetrahedralization.h
		src/tetwild/State.cpp
		src/tetwild/State.h
		src/tetwild/TetmeshElements.cpp
		src/tetwild/TetmeshElements.h
		src/tetwild/tetwild.cpp
		src/tetwild/VertexSmoother.cpp
		src/tetwild/VertexSmoother.h
		src/tetwild/geogram/mesh_AABB.cpp
		src/tetwild/geogram/mesh_AABB.h
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
		pymesh::pymesh
		spdlog::spdlog
	PRIVATE
		igl::cgal
		warnings::all
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
		igl::cgal
		warnings::all
)
target_include_directories(TetWild PRIVATE src)
igl_copy_cgal_dll(TetWild)

# Install
install(TARGETS TetWild RUNTIME DESTINATION bin)

################################################################################
# Folders for Visual Studio/XCode IDEs
################################################################################

# geogram
set_target_properties(geogram PROPERTIES FOLDER extern/geogram)
set_target_properties(geogram_third_party PROPERTIES FOLDER extern/geogram)
set_target_properties(uninstall PROPERTIES FOLDER extern/geogram)
# pymesh
set_target_properties(pymesh_tiny PROPERTIES FOLDER extern/pymesh)
