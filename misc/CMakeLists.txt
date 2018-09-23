add_executable(mmg_remesh mmg_remesh)
target_link_libraries(mmg_remesh
		tetwild::tetwild
		tetwild::internal
		CLI11::CLI11
		warnings::all
)
target_include_directories(mmg_remesh PRIVATE src)
igl_copy_cgal_dll(TetWild)
