cmake_minimum_required(VERSION 3.7.0...${CMAKE_VERSION})

if (NOT DEFINED BUILD_GPU)
SET (BUILD_GPU OFF)
endif (NOT DEFINED BUILD_GPU)

# --fmad=true and --fmad=false enables and disables the contraction respectively. 
 # This switch is supported only when the --gpu-architecture option is set with     
 # compute_20, sm_20, or higher. For other architecture classes, the contraction is     
  # always enabled. 
 # The --use_fast_math option implies --fmad=true, and enables the contraction.
if (BUILD_GPU)
set(CMAKE_CXX_STANDARD 11) #LIKE IN CUNSearch, this allows to include string (conflicts with CXX17)!!!
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fmad=false --prec-sqrt=true --ftz=false --prec-div=true")
set_source_files_properties(Domain_d.C Matrices.C Mechanical.C Solver.C Tensor.C main.C PROPERTIES LANGUAGE CUDA)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
ELSE ()
# if (MSVC)
    # # MSVC equivalent optimizations
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /std:c++17 /fp:precise")
# elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # # GCC, Clang, and MinGW
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -ffp-contract=off -fno-finite-math-only")
# endif()
ENDIF() 
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math=false")

# FROM --use_fast_math implies --ftz=true --prec-div=false --prec-sqrt=false --fmad=true.
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html


SET (BUILD_REMESH ON)
if (BUILD_OMEGA_H STREQUAL(ON))
message("SETTING OMEGA_H TO ON, FORCED REMESH TO ON")
else()
message("SETTING OMEGA_H TO OFF")
endif()


if (NOT DEFINED BUILD_OMEGA_H)
  set (BUILD_OMEGA_H OFF)
endif()

if (NOT DEFINED BUILD_MMG)
  set (BUILD_MMG OFF)
endif()

if (BUILD_MMG)
  add_compile_definitions(MMG_BUILD)
endif()

if (BUILD_GPU)
  add_compile_definitions(CUDA_BUILD)
else()
  # include_directories(${PROJECT_SOURCE_DIR}/lib/Math/include)
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(${PROJECT_SOURCE_DIR}/src/include)

MESSAGE ("CHECKING REMESH")
if (BUILD_REMESH)
    MESSAGE ("ADDING OMEGA_H")
  if (BUILD_OMEGA_H)
    set (Omega_h_USE_ZLIB OFF)
    MESSAGE ("ADDING OMEGA_H")
    if ((CMAKE_BUILD_TYPE STREQUAL "Debug"))
    message("OMEGA_H FORCED TO BUILD AS DEBUG")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    else()
    message("OMEGA_H FORCED TO BUILD AS RELEASE")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    
    if (BUILD_GPU)
    set (Omega_h_USE_CUDA ON)
    else()
      set (Omega_h_USE_OpenMP ON)
    #default is CPU
    endif() 
    add_subdirectory("lib/omega_h-9.34.13")
  endif()

  if (BUILD_MMG)
    MESSAGE ("ADDING MMG")
    add_subdirectory("lib/mmg-5.8.0")    
  endif()
#  ADD_CUSTOM_COMMAND(#TARGET _dnlPython 
#                  POST_BUILD
#                  COMMAND ${CMAKE_COMMAND} -E copy_if_different
#                  ${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src/Omega_h_config.h
#                  ${EXECUTABLE_OUTPUT_PATH})
  include_directories(${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h
  #add_subdirectory("examples/remesh")
endif (BUILD_REMESH)

#AFTER REMESH
add_subdirectory (src)

# if (BUILD_GPU)
# project (WeldFormFEM LANGUAGES CXX CUDA)
# else ()
# project (WeldFormFEM LANGUAGES CXX)
# endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(${PROJECT_SOURCE_DIR}/src/include)

