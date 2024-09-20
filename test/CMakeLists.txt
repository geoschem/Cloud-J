# Define a test that uses the path to the executable
add_test(NAME CompareOutput
  COMMAND ${CMAKE_COMMAND} -E env sh -c "$<TARGET_FILE:cloudj_standalone> > actual_output.txt && diff actual_output.txt ${CMAKE_SOURCE_DIR}/test/expected_output/reference_output.txt"
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
)
