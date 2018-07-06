################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
project(TetWild)

################################################################################
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################################

option(GTET_ISPC "Use ISPC" OFF)

# Build static binaries
#set(BUILD_SHARED_LIBS OFF)
set(BUILD_STATIC_RUNTIME OFF)
set(BUILD_DYNAMIC_RUNTIME ON)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#ADD_DEFINITIONS("-D_DLL")

################################################################################
# 3rd party libraries
################################################################################


# Geogram library
include(geogram)

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
find_package(LIBIGL REQUIRED QUIET)
find_package(Eigen REQUIRED)

################################################################################

SET(MAIN_FILE main.cpp)

set(SOURCE_FILES
		BSPElements.h TetmeshElements.h
		Preprocess.cpp Preprocess.h
		DelaunayTetrahedralization.cpp DelaunayTetrahedralization.h
		MeshConformer.cpp MeshConformer.h
		BSPSubdivision.cpp BSPSubdivision.h
		SimpleTetrahedralization.cpp SimpleTetrahedralization.h
		MeshRefinement.cpp MeshRefinement.h
		LocalOperations.cpp LocalOperations.h
		EdgeSplitter.cpp EdgeSplitter.h
		EdgeCollapser.cpp EdgeCollapser.h
		EdgeRemover.cpp EdgeRemover.h
		VertexSmoother.cpp VertexSmoother.h
		InoutFiltering.cpp InoutFiltering.h
		CLI11.hpp

		pymesh/MshLoader.cpp pymesh/MshLoader.h
		pymesh/MshSaver.cpp pymesh/MshSaver.h
		pymesh/Exception.h
		heads.cpp heads.h
		tetwild.h tetwild.cpp
		)


SET(OBJ_FILES "")
if(GTET_ISPC)
	add_custom_command(OUTPUT energy_ispc.o
			COMMAND ispc -g ${CMAKE_SOURCE_DIR}/ispc/energy.ispc -h ${CMAKE_SOURCE_DIR}/ispc/energy.h -o energy_ispc.o DEPENDS ${CMAKE_SOURCE_DIR}/ispc/energy.ispc)
	add_definitions(-DGTET_ISPC)
    LIST(APPEND OBJ_FILES energy_ispc.o)
endif(GTET_ISPC)
ADD_LIBRARY(lib_tetwild SHARED ${SOURCE_FILES} ${OBJ_FILES})
target_link_libraries(lib_tetwild geogram igl::core igl::cgal)
target_include_directories(lib_tetwild SYSTEM PUBLIC ${EIGEN_INCLUDE_DIRS})
set_target_properties(lib_tetwild PROPERTIES OUTPUT_NAME "TetWild")

add_executable(TetWild ${MAIN_FILE})
target_link_libraries(TetWild lib_tetwild)

#target_link_libraries(TetWild geogram igl::core igl::cgal igl::viewer)

INSTALL(TARGETS lib_tetwild TetWild
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)
INSTALL(FILES tetwild.h
    DESTINATION include)
