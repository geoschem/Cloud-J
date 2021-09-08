# Cloud-J/src/CMakeLists.txt

add_library(CloudJ_Standalone STATIC EXCLUDE_FROM_ALL
	CJ77.f90
)
target_link_libraries(CloudJ_Standalone
	PUBLIC CloudJ_Core
)
target_include_directories(CloudJ_Standalone
        INTERFACE ${CLOUDJ_BINARY_DIR}/mod
)
target_compile_options(CloudJ_Standalone
        PRIVATE
                ""
                $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)

#if("cloudj_standalone" IN_LIST CLOUDJ_EXE_TARGETS)
        add_executable(cloudj_standalone
                CJ77.f90
        )
        target_link_libraries(cloudj_standalone
                PUBLIC CloudJ_Standalone
        )
        set_target_properties(cloudj_standalone PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
        install(TARGETS cloudj_standalone
                RUNTIME DESTINATION ${RUNDIR}
        )
#endif()
