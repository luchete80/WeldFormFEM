CMAKE_MINIMUM_REQUIRED (VERSION 2.8)


if (BUILD_GPU)
project (MetFEMGPU LANGUAGES CXX CUDA)
project (WeldFormFEM LANGUAGES CXX CUDA)
project (testcuns LANGUAGES CXX CUDA)
else ()
project (MetFEMGPU LANGUAGES CXX)
project (WeldFormFEM LANGUAGES CXX)
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/cuda)
include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(../lib/LSDynaReader/src/)


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
set_source_files_properties(Domain_d.C Matrices.C Mechanical.C Solver.C Tensor.C main.C 
                            WeldFormFEM.C
                            PROPERTIES LANGUAGE CUDA)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
ELSE ()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
ENDIF() 
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math=false")

# FROM --use_fast_math implies --ftz=true --prec-div=false --prec-sqrt=false --fmad=true.
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html


# include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# set_source_files_properties(Domain_d.C PROPERTIES LANGUAGE CUDA)

add_executable(MetFEMGPU
                        main.C
                        Domain_d.C
                        Matrices.C
                        Mechanical.C
                        Matrices.C
                        Solver.C

              )
              
# add_executable(WeldFormFEM
                        # WeldFormFEM.C
                        # Domain_d.C
                        # Matrices.C
                        # Mechanical.C
                        # Matrices.C
                        # Solver.C
                        # ../lib/LSDynaReader/src/lsdynaReader.cpp

              # )
add_executable(1el_2d_sig_dev_pe  1el_2d_sig_dev_pe.cpp
                                  Domain_d.C
                                  Matrices.C
                                  Mechanical.C)

add_executable(1el_3d  test_1el_3D.cpp
                       Domain_d.C
                       Matrices.C
                       Mechanical.C)

#THIS SHOULD BE SET GREATER THAN 60 DUE TO CUDA ADDATOMIC
set_property(TARGET MetFEMGPU PROPERTY CUDA_ARCHITECTURES 60)

# set_property(TARGET WeldFormFEM PROPERTY CUDA_ARCHITECTURES 60)


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

if (BUILD_GPU)
TARGET_LINK_LIBRARIES (MetFEMGPU  cuda cudadevrt)
# TARGET_LINK_LIBRARIES (WeldFormFEM cuda cudadevrt)
else()
# TARGET_LINK_LIBRARIES (MetFEMGPU libstdc++)
endif()

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