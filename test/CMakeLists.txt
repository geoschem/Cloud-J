if(WIN32)
  add_test(NAME GenerateOutput
    COMMAND ${CMAKE_COMMAND} -E env cmd /c "$<TARGET_FILE:cloudj_standalone> > actual_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )

  add_test(NAME CompareOutput
    COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol actual_output.txt "${CMAKE_SOURCE_DIR}/test/expected_output/reference_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )
  set_tests_properties(CompareOutput PROPERTIES DEPENDS GenerateOutput)

  if(CLOUDJ_CPP_PORT)
    add_test(NAME GenerateOutputCpp
      COMMAND ${CMAKE_COMMAND} -E env cmd /c "$<TARGET_FILE:cloudj_standalone_cpp> --use-headers-only > cpp_actual_output.txt"
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
    )
    add_test(NAME CompareOutputCpp
      COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol cpp_actual_output.txt "${CMAKE_SOURCE_DIR}/test/expected_output/cpp_reference_output.txt"
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
    )
    set_tests_properties(CompareOutputCpp PROPERTIES DEPENDS GenerateOutputCpp)
  endif()
else()
  add_test(NAME CompareOutput
    COMMAND ${CMAKE_COMMAND} -E env sh -c "ulimit -s hard 2>/dev/null; $<TARGET_FILE:cloudj_standalone> > actual_output.txt; ${CMAKE_COMMAND} -E compare_files --ignore-eol actual_output.txt ${CMAKE_SOURCE_DIR}/test/expected_output/reference_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )

  if(CLOUDJ_CPP_PORT)
    add_test(NAME CompareOutputCpp
      COMMAND ${CMAKE_COMMAND} -E env sh -c "$<TARGET_FILE:cloudj_standalone_cpp> --use-headers-only > cpp_actual_output.txt; ${CMAKE_COMMAND} -E compare_files --ignore-eol cpp_actual_output.txt ${CMAKE_SOURCE_DIR}/test/expected_output/cpp_reference_output.txt"
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
    )
  endif()
endif()

if(CLOUDJ_CPP_PORT)
  add_test(NAME FuzzPropertyParity
      COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/tools/run_property_fuzzer.py --count 10000
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )
endif()

if(CLOUDJ_CPP_PORT)
  add_executable(test_cpp_math
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_cpp_math.cpp
  )
  target_compile_definitions(test_cpp_math PRIVATE CLOUDJ_GPU_MODE)
  set_target_properties(test_cpp_math PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  )
  if(CLOUDJ_USE_KOKKOS)
    target_link_libraries(test_cpp_math PUBLIC Kokkos::kokkos Kokkos::kokkoskernels)
  endif()

  add_executable(test_fortran_math
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/test_fortran_math.F90
  )
  target_link_libraries(test_fortran_math PRIVATE CloudJ_Core)
  set_target_properties(test_fortran_math PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  )

  add_executable(test_osa ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_osa.cpp)
  set_target_properties(test_osa PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestOSA COMMAND test_osa)

  add_executable(test_safe_div ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_safe_div.cpp)
  set_target_properties(test_safe_div PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestSafeDiv COMMAND test_safe_div)

  add_executable(test_init ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_init.cpp)
  target_compile_definitions(test_init PRIVATE MODEL_STANDALONE)
  set_target_properties(test_init PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestInit COMMAND test_init)

  add_executable(test_extral1 ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_extral1.cpp)
  target_compile_definitions(test_extral1 PRIVATE MODEL_STANDALONE)
  set_target_properties(test_extral1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestExtral1 COMMAND test_extral1)

  add_executable(test_sphere ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_sphere.cpp)
  target_compile_definitions(test_sphere PRIVATE MODEL_STANDALONE)
  set_target_properties(test_sphere PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestSphere COMMAND test_sphere)

  add_executable(test_optics ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_optics.cpp)
  target_compile_definitions(test_optics PRIVATE MODEL_STANDALONE)
  set_target_properties(test_optics PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestOptics COMMAND test_optics)

  add_executable(test_photo_jx ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_photo_jx.cpp)
  target_compile_definitions(test_photo_jx PRIVATE MODEL_STANDALONE)
  set_target_properties(test_photo_jx PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestPhotoJX COMMAND test_photo_jx)

  add_executable(test_heapsort ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_heapsort.cpp)
  target_compile_definitions(test_heapsort PRIVATE MODEL_STANDALONE)
  set_target_properties(test_heapsort PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestHeapsort COMMAND test_heapsort)

  add_executable(test_cloud_jx ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_cloud_jx.cpp)
  target_compile_definitions(test_cloud_jx PRIVATE MODEL_STANDALONE)
  set_target_properties(test_cloud_jx PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_test(NAME TestCloudJX COMMAND test_cloud_jx)

  add_executable(test_library_api ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_library_api.cpp)
  target_compile_definitions(test_library_api PRIVATE MODEL_STANDALONE
    CLOUDJ_TEST_TABLES_DIR="${CMAKE_BINARY_DIR}/bin/tables"
    CLOUDJ_TEST_GOLDEN_DIR="${CMAKE_SOURCE_DIR}/test/expected_output"
    CLOUDJ_TEST_FORTRAN7="${CMAKE_SOURCE_DIR}/bin/fort.7")
  set_target_properties(test_library_api PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  add_dependencies(test_library_api copy_table_data)
  add_test(NAME TestLibraryAPI COMMAND test_library_api)
endif()



