# Cloud-J/src/CMakeLists.txt

#if("cloudj_standalone" IN_LIST CLOUDJ_EXE_TARGETS)
        add_executable(cloudj_standalone
                cldj_standalone.F90
        )
        target_link_libraries(cloudj_standalone
                PUBLIC CloudJ_Core
        )
        set_target_properties(cloudj_standalone PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
        install(TARGETS cloudj_standalone
                RUNTIME DESTINATION ${RUNDIR}
        )
#endif()
