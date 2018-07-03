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
option(BUILD_WITH_DEBUG             "Debug mode"         OFF)
option(BUILD_AS_LIB                 "Build as library"   OFF)

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

################################################################################
# TetWild
################################################################################
include_directories(include)

set(SOURCE_FILES
	src/Preprocess.cpp
	src/DelaunayTetrahedralization.cpp
	src/MeshConformer.cpp
	src/BSPSubdivision.cpp
	src/SimpleTetrahedralization.cpp
	src/MeshRefinement.cpp
	src/LocalOperations.cpp
	src/EdgeSplitter.cpp
	src/EdgeCollapser.cpp
	src/EdgeRemover.cpp
	src/VertexSmoother.cpp
	src/InoutFiltering.cpp
	src/heads.cpp
	src/tetwild.cpp
	)

# Build with PyMesh for debug mode or when building an executable
if(BUILD_WITH_DEBUG OR NOT BUILD_AS_LIB)
	message(STATUS "Building with PyMesh.")
	set(SOURCE_FILES ${SOURCE_FILES}
		pymesh/MshLoader.cpp
		pymesh/MshSaver.cpp
		)
	include_directories(.)
    add_definitions(-DUSE_PYMESH)
endif()

# Define debug build
if(BUILD_WITH_DEBUG)
    add_definitions(-DDEBUG_MODE -DUSE_PYMESH)
endif()

# ispc
SET(OBJ_FILES "")
if(GTET_ISPC)
	add_custom_command(OUTPUT energy_ispc.o
			COMMAND ispc ${CMAKE_SOURCE_DIR}/ispc/energy.ispc -h ${CMAKE_SOURCE_DIR}/ispc/energy.h -o energy_ispc.o DEPENDS ${CMAKE_SOURCE_DIR}/ispc/energy.ispc)
	add_definitions(-DGTET_ISPC)
    LIST(APPEND OBJ_FILES energy_ispc.o)
endif(GTET_ISPC)

# Build shared library otherwise build a static library for an executable
if(BUILD_AS_LIB)
    message(STATUS "Building a shared library.")
    set(BUILD_SHARED_LIBS ON)
else()
    message(STATUS "Building a static library.")
    set(BUILD_SHARED_LIBS OFF)
endif()

add_library(libTetWild ${SOURCE_FILES} ${OBJ_FILES})
target_link_libraries(libTetWild geogram igl::core igl::cgal)
target_include_directories(libTetWild SYSTEM PUBLIC ${EIGEN_INCLUDE_DIRS})
if(BUILD_AS_LIB)
    set_target_properties(libTetWild PROPERTIES OUTPUT_NAME "tetwild")
endif()
if(WIN32)
    set_target_properties(libTetWild PROPERTIES COMPILE_FLAGS "-DWIN32 -DTETWILD_EXPORTS")
endif()

# Building exectuable
if(NOT BUILD_AS_LIB)
    message(STATUS "Building an executable.")
    set(MAIN_FILE src/main.cpp)
    add_executable(TetWild ${MAIN_FILE})
    target_link_libraries(TetWild libTetWild)
endif()

# Install
if(BUILD_AS_LIB)
    INSTALL(TARGETS libTetWild
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)
    INSTALL(FILES include/tetwild.h DESTINATION include/tetwild)
else()
    INSTALL(TARGETS TetWild RUNTIME DESTINATION bin)
endif()
