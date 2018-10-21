# - Try to find the LIBIGL library
# Once done this will define
#
#  LIBIGL_FOUND - system has LIBIGL
#  LIBIGL_INCLUDE_DIR - **the** LIBIGL include directory
if(LIBIGL_FOUND OR TARGET igl::core)
    return()
endif()

find_path(LIBIGL_INCLUDE_DIR igl/readOBJ.h
    HINTS
        # ENV LIBIGL
        # ENV LIBIGLROOT
        # ENV LIBIGL_ROOT
        # ENV LIBIGL_DIR
    PATHS
        ${TETWILD_EXTERNAL}/libigl
        # /usr
        # /usr/local
        # /usr/local/igl/libigl
    PATH_SUFFIXES include
    NO_CMAKE_SYSTEM_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBIGL
    "\nlibigl not found --- You can download it using:\n\tgit clone --recursive https://github.com/libigl/libigl.git ${CMAKE_SOURCE_DIR}/../libigl"
    LIBIGL_INCLUDE_DIR)
mark_as_advanced(LIBIGL_INCLUDE_DIR)

#list(APPEND CMAKE_MODULE_PATH "${LIBIGL_INCLUDE_DIR}/../shared/cmake")
list(APPEND CMAKE_MODULE_PATH "${LIBIGL_INCLUDE_DIR}/../cmake")
include(libigl)
