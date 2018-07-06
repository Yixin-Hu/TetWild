cmake_minimum_required(VERSION 3.10)

if(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267")
endif()

### Sources ###

set(sources
    heads.cpp
    heads.h
    EdgeSplitter.cpp
    EdgeSplitter.h
    EdgeCollapser.cpp
    EdgeCollapser.h
    EdgeRemover.cpp
    EdgeRemover.h
    tetwild.cpp
    tetwild.h
    LocalOperations.cpp
    LocalOperations.h
    Preprocess.cpp
    Preprocess.h
    MeshConformer.cpp
    MeshConformer.h
    DelaunayTetrahedralization.cpp
    DelaunayTetrahedralization.h
    BSPSubdivision.cpp
    BSPSubdivision.h
    SimpleTetrahedralization.cpp
    SimpleTetrahedralization.h
    MeshRefinement.cpp
    MeshRefinement.h
    VertexSmoother.cpp
    VertexSmoother.h
    InoutFiltering.cpp
    InoutFiltering.h
    pymesh/MshSaver.cpp
    pymesh/MshSaver.h
    pymesh/MshLoader.cpp
    pymesh/MshLoader.h
)

# This corrects cmake's default Visual Studio's ugly folder organization
source_group("" FILES ${sources})

add_library(tetwild ${sources})

target_include_directories(tetwild PUBLIC .)

target_link_libraries(tetwild PUBLIC
  glm geogram libigl2 cgal2)
