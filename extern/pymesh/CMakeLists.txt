add_library(pymesh_tiny
		MshLoader.h
		MshLoader.cpp
		MshSaver.h
		MshSaver.cpp
		Exception.h
)
target_include_directories(pymesh_tiny PUBLIC ..)

# Link against igl_common to draw Eigen and C++11
target_link_libraries(pymesh_tiny PUBLIC igl_common)

add_library(pymesh::pymesh ALIAS pymesh_tiny)
