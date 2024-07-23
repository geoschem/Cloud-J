# Cloud-J/src/Core/CMakeLists.txt

add_library(CloudJ_Core STATIC EXCLUDE_FROM_ALL
        cldj_cmn_mod.F90
        cldj_error_mod.F90
        cldj_init_mod.F90
        cldj_sub_mod.F90
        cldj_fjx_sub_mod.F90
        cldj_fjx_osa_mod.F90
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
