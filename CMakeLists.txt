cmake_minimum_required(VERSION 3.2)


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
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-fast-math -std=c++17 -O3")
set_source_files_properties(Domain_d.C Matrices.C Mechanical.C Solver.C Tensor.C main.C PROPERTIES LANGUAGE CXX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
ENDIF() 
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math=false")

# FROM --use_fast_math implies --ftz=true --prec-div=false --prec-sqrt=false --fmad=true.
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html

if (NOT DEFINED BUILD_GPU)
SET (BUILD_GPU OFF)
SET (BUILD_LAGRANGIAN ON)
SET (BUID_OMEGA_H ON) #SET OFF IFALREADY BUILT ANDYOURE TESTING
endif (NOT DEFINED BUILD_GPU)

if (BUILD_GPU)
  add_compile_definitions(CUDA_BUILD)
else()
  # include_directories(${PROJECT_SOURCE_DIR}/lib/Math/include)
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(${PROJECT_SOURCE_DIR}/src/include)



add_subdirectory (src)

if (BUILD_LAGRANGIAN)
  if (BUID_OMEGA_H)
    add_subdirectory("lib/omega_h-9.34.13")
    set (Omega_h_USE_ZLIB OFF)
    if (BUILD_GPU)
    set (Omega_h_USE_CUDA ON)
    else()
      set (Omega_h_USE_OpenMP ON)
    #default is CPU
    endif() 
  endif()
  add_subdirectory("examples/remesh")
endif (BUILD_LAGRANGIAN)

if (BUILD_GPU)
project (WeldFormFEM LANGUAGES CXX CUDA)
else ()
project (WeldFormFEM LANGUAGES CXX)
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(${PROJECT_SOURCE_DIR}/src/include)



link_directories(${CUNSEARCH_LIB_DIR})
link_directories(${cuNSearch_BINARY_DIR})
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
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --fmad=false -std=c++17 -O3 -fno-finite-math-only")
ENDIF() 
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math=false")

# FROM --use_fast_math implies --ftz=true --prec-div=false --prec-sqrt=false --fmad=true.
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html



              
# add_executable(WeldFormFEM
                        # WeldForm.C
                        # Domain_d.C
                        # Matrices.C
                        # Mechanical.C
                        # Matrices.C
                        # Solver.C

              # )



