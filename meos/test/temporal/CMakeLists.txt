if(NOT MEOS_EXPERIMENTAL_ANALYTICS)
  return()
endif()

add_executable(test_tekf_gps
  test_tekf_gps.c
)

target_include_directories(test_tekf_gps PRIVATE
  ${CMAKE_SOURCE_DIR}/meos/include
  ${CMAKE_BINARY_DIR}
)

target_link_libraries(test_tekf_gps PRIVATE meos)
target_compile_definitions(test_tekf_gps PRIVATE MEOS_EXPERIMENTAL_ANALYTICS=1)

add_test(NAME meos_test_tekf_gps COMMAND test_tekf_gps)
