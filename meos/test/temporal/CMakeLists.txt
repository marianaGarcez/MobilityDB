if(NOT MEOS_EXPERIMENTAL_ANALYTICS)
  return()
endif()

add_executable(test_tkalman
  test_tkalman.c
)

target_include_directories(test_tkalman PRIVATE
  ${CMAKE_SOURCE_DIR}/meos/include
  ${CMAKE_BINARY_DIR}
)

target_link_libraries(test_tkalman PRIVATE meos)

# Ensure the prototype is visible in the test TU
target_compile_definitions(test_tkalman PRIVATE MEOS_EXPERIMENTAL_ANALYTICS=1)

add_test(NAME meos_test_tkalman COMMAND test_tkalman)

add_executable(bench_tkalman
  bench_tkalman.c
)

target_include_directories(bench_tkalman PRIVATE
  ${CMAKE_SOURCE_DIR}/meos/include
  ${CMAKE_BINARY_DIR}
)

target_link_libraries(bench_tkalman PRIVATE meos)
target_compile_definitions(bench_tkalman PRIVATE MEOS_EXPERIMENTAL_ANALYTICS=1)

add_test(NAME meos_bench_tkalman COMMAND bench_tkalman)
