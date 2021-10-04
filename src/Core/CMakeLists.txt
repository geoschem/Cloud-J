# Cloud-J/src/Core/CMakeLists.txt

add_library(CloudJ_Core STATIC EXCLUDE_FROM_ALL
        cld_sub_mod.f90
        fjx_error_mod.f90
        fjx_cmn_mod.f90
        fjx_sub_mod.f90
        fjx_init_mod.f90
        osa_sub_mod.f90
)
target_link_libraries(CloudJ_Core
	PUBLIC CloudJBuildProperties
)
target_include_directories(CloudJ_Core
	INTERFACE ${CLOUDJ_BINARY_DIR}/mod
)
target_compile_options(CloudJ_Core
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)
