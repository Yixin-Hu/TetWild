################################################################################
# CMake download helpers
################################################################################

# download external dependencies
include(TetWildDownloadExternal)

################################################################################
# Required dependencies
################################################################################

# Boost
if(TETWILD_WITH_HUNTER)
	hunter_add_package(Boost COMPONENTS thread system)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
	tetwild_download_spdlog()
	add_subdirectory(${TETWILD_EXTERNAL}/spdlog)
endif()

# libigl
if(NOT TARGET igl::core)
	tetwild_download_libigl()
	find_package(LIBIGL REQUIRED)
endif()

# geogram
if(NOT TARGET geogram)
	tetwild_download_geogram()
	include(geogram)
endif()

# pymesh loaders
add_subdirectory(${TETWILD_EXTERNAL}/pymesh)

# CL11
if(NOT TARGET CLI11::CLI11)
	tetwild_download_cli11()
	add_subdirectory(${TETWILD_EXTERNAL}/cli11)
endif()
