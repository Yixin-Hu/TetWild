################################################################################
# CMake download helpers
################################################################################

# download external dependencies
include(TetWildDownloadExternal)

################################################################################
# Required dependencies
################################################################################

# mmg
if(NOT TARGET mmg)
	tetwild_download_mmg()
	option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
	option(BUILD_TESTING "Enable/Disable continuous integration" OFF)
	set(OLD_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS} CACHE BOOL "" FORCE)
	set(OLD_BUILD_TESTING ${BUILD_TESTING} CACHE BOOL "" FORCE)
	set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
	set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
	add_subdirectory(${TETWILD_EXTERNAL}/mmg)
	set(BUILD_SHARED_LIBS ${OLD_BUILD_SHARED_LIBS} CACHE BOOL "" FORCE)
	set(BUILD_TESTING ${OLD_BUILD_TESTING} CACHE BOOL "" FORCE)
	add_library(mmg::mmg ALIAS libmmg_so)
	add_library(mmg::mmgs ALIAS libmmgs_so)
	add_library(mmg::mmg2d ALIAS libmmg2d_so)
	add_library(mmg::mmg3d ALIAS libmmg3d_so)

	foreach(target_name IN ITEMS libmmg_so libmmgs_so libmmg2d_so libmmg3d_so)
		set_target_properties(${target_name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
		set_target_properties(${target_name} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
		set_target_properties(${target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
	endforeach()

endif()
