add_library(pymesh_tiny
		MshLoader.h
		MshLoader.cpp
		MshSaver.h
		MshSaver.cpp
		Exception.h
)
target_include_directories(pymesh_tiny PUBLIC ..)

if(TARGET Eigen3::Eigen)
	# If an imported target already exists, use it
	target_link_libraries(pymesh_tiny PUBLIC Eigen3::Eigen)
else()
	# Otherwise use the version bundled with libigl
	target_include_directories(pymesh_tiny SYSTEM PUBLIC ${LIBIGL_EXTERNAL}/eigen)
endif()

add_library(pymesh::pymesh ALIAS pymesh_tiny)
