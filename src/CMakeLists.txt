cmake_minimum_required(VERSION 3.2)

if (NOT DEFINED BUILD_GPU)
SET(BUILD_GPU OFF)
endif (NOT DEFINED BUILD_GPU)

if (BUILD_GPU)
  add_compile_definitions(CUDA_BUILD)
endif (BUILD_GPU)

PROJECT (MetFEMGPU)

if (BUILD_GPU)
project (MetFEMGPU LANGUAGES CXX CUDA)
project (testcuns LANGUAGES CXX CUDA)
else ()
project (MetFEMGPU LANGUAGES CXX)
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/cuda)
include_directories(${PROJECT_SOURCE_DIR}/include)

link_directories(${CUNSEARCH_LIB_DIR})
link_directories(${cuNSearch_BINARY_DIR})

set(CMAKE_CXX_STANDARD 11) #LIKE IN CUNSearch, this allows to include string (conflicts with CXX17)!!!
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --fmad=true and --fmad=false enables and disables the contraction respectively. 
 # This switch is supported only when the --gpu-architecture option is set with     
 # compute_20, sm_20, or higher. For other architecture classes, the contraction is     
  # always enabled. 
 # The --use_fast_math option implies --fmad=true, and enables the contraction.
if (BUILD_GPU)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fmad=false --prec-sqrt=true --ftz=false --prec-div=true")
set_source_files_properties(Domain_d.C PROPERTIES LANGUAGE CUDA)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
ELSE ()

ENDIF() 
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math=false")

# FROM --use_fast_math implies --ftz=true --prec-div=false --prec-sqrt=false --fmad=true.
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html


# include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# set_source_files_properties(Domain_d.C PROPERTIES LANGUAGE CUDA)

add_executable(MetFEMGPU
                          main.cu
													Domain_d.C
                          cuda/Matrices.cu
                          cuda/Mechanical.cu
                          cuda/cudautils.cu
													cuda/Solver.cu

              )

#THIS SHOULD BE SET GREATER THAN 60 DUE TO CUDA ADDATOMIC
set_property(TARGET MetFEMGPU PROPERTY CUDA_ARCHITECTURES 60)



# add_executable(main_thermal 
													# main_thermal.cu
													# Domain.cpp 
													# #NbSearch.cpp 
													# Output.cpp
													# Vector.cpp
                          # cuda/thermal.cu
													# #TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
													# cuda/Domain_d.cu
													# #cuda/Thermal.cu
													# #cuda/Interaction_global.cu	#TEMPORARY, ALL IS EXTERNAL
													# cuNSearch/src/PointSet.cu cuNSearch/src/PointSetImplementation.cu
													# cuNSearch/src/cuNSearch.cu cuNSearch/src/cuNSearchDeviceData.cu cuNSearch/src/cuNSearchKernels.cu
													# cuNSearch/Utils/cuda_helper.cu cuNSearch/Utils/Timing.cpp
													# cuNSearch/Utils/IDFactory.cpp													
													# )
                          

# add_executable(main_pairs 
													# #main_thermal.cu
													# main_thermal_pairs.cu
													# Domain.cpp 
													# #NbSearch.cpp 
													# Output.cpp
													# Vector.cpp
													# #TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
													# cuda/Domain_d.cu
													# cuda/Thermal.cu
													# #cuda/Interaction_global.cu	#TEMPORARY, ALL IS EXTERNAL
													# cuNSearch/src/PointSet.cu
													# cuNSearch/src/PointSetImplementation.cu
													# cuNSearch/src/cuNSearch.cu
													# cuNSearch/src/cuNSearchDeviceData.cu
													# cuNSearch/src/cuNSearchKernels.cu
													# cuNSearch/Utils/cuda_helper.cu
													# cuNSearch/Utils/Timing.cpp
													# cuNSearch/Utils/IDFactory.cpp													
													# ) 


#tensor and all functions defined in separate modulus
# CASE NO INLINE VERSION
#set_property(TARGET gpuweldform PROPERTY CUDA_SEPARABLE_COMPILATION ON)

#add_subdirectory (cuda)

TARGET_LINK_LIBRARIES (MetFEMGPU cuda cudadevrt)

#FOREACH(var ${EXES})
#    ADD_EXECUTABLE        (${var} "${var}.cpp" #
#
#	#"${CMAKE_SOURCE_DIR}/Source/CompactNSearch/CompactNSearch.cpp"
#	)
#    TARGET_LINK_LIBRARIES (${var} ${LIBS} )
#    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
#ENDFOREACH(var)

# if (NOT BUILD_GPU)
# else()
# TARGET_LINK_LIBRARIES (${var} cuda)
# endif()
    # ADD_EXECUTABLE        (${EXES} "Plate-Yield.cpp")
    # TARGET_LINK_LIBRARIES (${EXES} ${LIBS} blitz)
    # SET_TARGET_PROPERTIES (${EXES} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
# else()
# add_executable(MetFEMGPU
              # main.cpp)

# endif ()