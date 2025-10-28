add_library(analytics OBJECT
  tkalman.c
)

target_include_directories(analytics PRIVATE
  ${CMAKE_SOURCE_DIR}/meos/include
  ${CMAKE_BINARY_DIR}
)
