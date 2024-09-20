# Define a test that uses the path to the executable
if(WIN32)
  add_test(NAME CompareOutput
    COMMAND ${CMAKE_COMMAND} -E env cmd /c "$<TARGET_FILE:cloudj_standalone> > actual_output.txt && ${CMAKE_COMMAND} -E compare_files actual_output.txt ${CMAKE_SOURCE_DIR}/test/expected_output/reference_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )
else()
  add_test(NAME CompareOutput
    COMMAND ${CMAKE_COMMAND} -E env sh -c "$<TARGET_FILE:cloudj_standalone> > actual_output.txt && ${CMAKE_COMMAND} -E compare_files actual_output.txt ${CMAKE_SOURCE_DIR}/test/expected_output/reference_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )
endif()