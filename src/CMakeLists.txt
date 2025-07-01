add_subdirectory(common)
include_directories(../include/common)
if (BUILD_EXPLICIT)
  add_subdirectory(explicit)
endif()

if (BUILD_IMPLICIT)
  add_subdirectory(implicit)
endif()
