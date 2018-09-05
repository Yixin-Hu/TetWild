set(tetwild_ispc__internal_dir ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")

function(ispc_add_energy target_name)
	# Compilation flags for ISPC
	if(CMAKE_BUILD_TYPE MATCHES Release)
		set(TETWILD_ISPC_FLAGS "")
	else()
		set(TETWILD_ISPC_FLAGS "-g")
	endif()

	add_custom_command(
		COMMAND
			ispc --pic --target=host ${TETWILD_ISPC_FLAGS}
				${tetwild_ispc__internal_dir}/energy.ispc
				-h ${tetwild_ispc__internal_dir}/energy.h
				-o ${CMAKE_CURRENT_BINARY_DIR}/energy_ispc.o
		DEPENDS
			${tetwild_ispc__internal_dir}/energy.ispc
		OUTPUT
			${tetwild_ispc__internal_dir}/energy.h
			${CMAKE_CURRENT_BINARY_DIR}/energy_ispc.o
	)

	target_sources(${target_name} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/energy_ispc.o)
	target_compile_definitions(${target_name} PUBLIC -DTETWILD_WITH_ISPC)
	target_include_directories(${target_name} PUBLIC ${tetwild_ispc__internal_dir}/..)
endfunction()
