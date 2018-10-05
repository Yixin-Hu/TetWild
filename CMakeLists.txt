cmake_minimum_required(VERSION 3.3)
project(TetWild)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4005 /wd4244 /wd4028 /wd4267 /wd4101 /wd4334 /wd4018 /wd4717 /wd4305")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4005 /wd4244 /wd4028 /wd4267 /wd4101 /wd4334 /wd4018 /wd4717 /wd4305")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(TETWILD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/extern")

include(TetWildDependencies)

################################################################################
# TetWild
################################################################################

add_library(tetwild STATIC
		extern/pymesh/MshLoader.h
		extern/pymesh/MshLoader.cpp
		extern/pymesh/MshSaver.h
		extern/pymesh/MshSaver.cpp
		extern/pymesh/Exception.h
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
		src/tetwild/Quality.cpp
		src/tetwild/Quality.h
		src/tetwild/SimpleTetrahedralization.cpp
		src/tetwild/SimpleTetrahedralization.h
		src/tetwild/State.cpp
		src/tetwild/State.h
		src/tetwild/TetmeshElements.cpp
		src/tetwild/TetmeshElements.h
		src/tetwild/tetwild.cpp
		src/tetwild/Utils.cpp
		src/tetwild/Utils.h
		src/tetwild/VertexSmoother.cpp
		src/tetwild/VertexSmoother.h
		src/tetwild/geogram/MeshAABB.cpp
		src/tetwild/geogram/MeshAABB.h
		src/tetwild/geogram/Utils.cpp
		src/tetwild/geogram/Utils.h
		src/tetwild/mmg/Remeshing.cpp
		src/tetwild/mmg/Remeshing.h
		src/tetwild/mmg/internal/Conversion.cpp
		src/tetwild/mmg/internal/Conversion.h
)
target_include_directories(tetwild
	PUBLIC
		extern
		src
		include
)
target_link_libraries(tetwild
	PRIVATE
		libigl2
		cgal2
		geogram
	PUBLIC
		spdlog
		mmg::mmg
)

# ispc
if(TETWILD_WITH_ISPC)
	message(STATUS "Compiling energy with ISPC")
	add_subdirectory(src/ispc)
	ispc_add_energy(tetwild)
endif()
