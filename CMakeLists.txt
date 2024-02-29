cmake_minimum_required(VERSION 3.2)

if (NOT DEFINED BUILD_GPU)
SET(BUILD_GPU OFF)
endif (NOT DEFINED BUILD_GPU)

if (BUILD_GPU)
  add_compile_definitions(CUDA_BUILD)
else()
  # include_directories(${PROJECT_SOURCE_DIR}/lib/Math/include)
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(${PROJECT_SOURCE_DIR}/src/include)



add_subdirectory (${PROJECT_SOURCE_DIR}/src)

