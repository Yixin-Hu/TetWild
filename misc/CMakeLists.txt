add_executable(levelset levelset.cpp)
target_link_libraries(levelset
		tetwild::tetwild
		tetwild::internal
		CLI11::CLI11
		warnings::all
)
target_include_directories(levelset PRIVATE src)
igl_copy_cgal_dll(levelset)

if(TETWILD_WITH_SANITIZERS)
	add_sanitizers(levelset)
endif()
