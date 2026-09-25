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
        if(WIN32)
                target_link_libraries(cloudj_standalone
                        PRIVATE "-Wl,--stack,10485760"
                )
        endif()
#endif()

if(CLOUDJ_CPP_PORT)
        add_executable(cloudj_standalone_cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/standalone/cloudj_standalone.cpp
        )
        set_target_properties(cloudj_standalone_cpp PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
        )

        option(CLOUDJ_GPU_MODE_BENCH "Compile C++ standalone in GPU Hermite spline mode" OFF)
        if(CLOUDJ_GPU_MODE_BENCH)
                target_compile_definitions(cloudj_standalone_cpp PRIVATE CLOUDJ_GPU_MODE)
        endif()

        if(WIN32)
                target_link_libraries(cloudj_standalone_cpp
                        PRIVATE "-Wl,--stack,10485760"
                )
        endif()

        if(CLOUDJ_USE_KOKKOS)
                target_link_libraries(cloudj_standalone_cpp PUBLIC Kokkos::kokkos Kokkos::kokkoskernels)
        endif()

        if(CLOUDJ_USE_OPENMP)
                target_link_libraries(cloudj_standalone_cpp PUBLIC OpenMP::OpenMP_CXX)
        endif()

        option(CLOUDJ_CXX_MARCH_NATIVE "Enable -march=native for C++ targets (opt-in; not portable)" OFF)
        if(CLOUDJ_CXX_MARCH_NATIVE)
                target_compile_options(cloudj_standalone_cpp PRIVATE -march=native)
        endif()

        # Match the Fortran Release optimisation flags (CloudJ_Fortran_FLAGS_-
        # RELEASE_GNU = -O3 -funroll-loops): the hot solver loops are fixed-size
        # 4x4 / 8-wide and benefit from unrolling the same way the Fortran does.
        target_compile_options(cloudj_standalone_cpp PRIVATE
                $<$<CONFIG:Release>:-funroll-loops>
                $<$<CONFIG:RelWithDebInfo>:-funroll-loops>
        )
endif()

