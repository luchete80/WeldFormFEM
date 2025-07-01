CMAKE_MINIMUM_REQUIRED (VERSION 3.5)

set(PROJECT_VERSION "0.0.1")  # or your current version
message(STATUS "Version: ${PROJECT_VERSION}")

set (BUILD_EXAMPLES OFF)


if (BUILD_GPU)
project (WeldFormFEM LANGUAGES CXX CUDA )
else ()
project (WeldFormFEM LANGUAGES CXX )
endif (BUILD_GPU)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/cuda)
include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(../lib/LSDynaReader/src/)

include_directories(../lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h
set (MMG_ROOT "${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/")
link_directories("${MMG_ROOT}/lib")
include_directories("${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/include") #FOR: Omega_h_config.h

set(OMEGA_H_DIR "${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src")  # Adjust this path
link_directories(${OMEGA_H_DIR})

set(CMAKE_PREFIX_PATH "/path/to/omega_h")
include_directories(${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h


# --fmad=true and --fmad=false enables and disables the contraction respectively. 
 # This switch is supported only when the --gpu-architecture option is set with     
 # compute_20, sm_20, or higher. For other architecture classes, the contraction is     
  # always enabled. 
 # The --use_fast_math option implies --fmad=true, and enables the contraction.
if (BUILD_GPU)
set(CMAKE_CXX_STANDARD 11) #LIKE IN CUNSearch, this allows to include string (conflicts with CXX17)!!!
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fmad=false --prec-sqrt=true --ftz=false --prec-div=true")
set_source_files_properties(Contact.C Domain_d.C Matrices.C Mechanical.C input.C Mesh.C Solver.C Tensor.C main.C  main_1_elem_3d.C
                            WeldFormFEM.C VTKWriter.C 
                            include/*.h src/include/defs.h src/include/utils.h
                            src/include/*.h
                            ${PROJECT_SOURCE_DIR}/src/include/*.cuh
                            main_4_elem_3d.C
                            main_1_elem_3d.C
                            PROPERTIES LANGUAGE CUDA)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
ELSE ()
  INCLUDE      (FindOpenMP)
  message ("BUILD GPU set to OFF")
  if (OPENMP_FOUND)
      message ("OPENMP FOUND!!!")
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else()
  message ("OPENMP NOT FOUND")
  endif()
set_source_files_properties(Domain_d.C Matrices.C Mechanical.C Solver.C Tensor.C main.C  main_1_elem_3d.C
                            ReMesher.C ReMesher_mmg.C
                            WeldFormFEM.C include/*.h src/include/defs.h src/include/utils.h
                            ${PROJECT_SOURCE_DIR}/src/include/*.cuh
                            main_4_elem_3d.C
                            PROPERTIES LANGUAGE CXX)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
ENDIF() 
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math=false")

# FROM --use_fast_math implies --ftz=true --prec-div=false --prec-sqrt=false --fmad=true.
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html


# include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# set_source_files_properties(Domain_d.C PROPERTIES LANGUAGE CUDA)

# add_executable(MetFEMGPU
                        # main.C
                        # Domain_d.C
                        # Matrices.C
                        # Mechanical.C
                        # Matrices.C
                        # Solver.C
                        # VTKWriter.C

              # )
 #add_executable(1elem_2d
#                         main_1_elem_2d.C
#                         Contact.C
#                         Domain_d.C
#                         Matrices.C
#                         Mechanical.C
#                         Matrices.C
#                         Solver.C

#                         VTKWriter.C

#               )


set (SOURCES                          
                        Contact.C
                         Domain_d.C
                         Matrices.C
                         Mechanical.C
                         Mesh.C
                         Matrices.C
                         Solver_explicit.C
                         Thermal.C
                         input.C
                         ../lib/LSDynaReader/src/lsdynaReader.cpp
                         VTKWriter.C
                   )

if (BUILD_REMESH)
 set (SOURCES "${SOURCES}" ReMesher.C
                           ReMesher_mmg.C)
endif()

if (BUILD_EXAMPLES)

 add_executable(1elem
                         main_1_elem_3d.C
                         "${SOURCES}"

               )

 add_executable(4elem
                         main_4_elem_3d.C
                         "${SOURCES}"

               )

endif()
                   
add_executable(WeldFormFEM WeldFormFEM.C
                         "${SOURCES}"
               )

if (BUILD_OMEGA_H)
  #add_executable (test_omega_h test_omega_h.cpp)
endif()


if (BUILD_GPU)

set_property(TARGET WeldFormFEM PROPERTY CUDA_ARCHITECTURES 75)


if (BUILD_EXAMPLES)
#TARGET_LINK_LIBRARIES (1elem  cuda cudadevrt)
TARGET_LINK_LIBRARIES (1elem   cuda cudart cudadevrt)
TARGET_LINK_LIBRARIES (4elem  cuda cudart cudadevrt)
endif()

TARGET_LINK_LIBRARIES (WeldFormFEM cuda cudart cudadevrt)
else() #BUILD CPU
  # TARGET_LINK_LIBRARIES (MetFEMGPU libstdc++)
  if (BUILD_REMESH)
  if (BUILD_EXAMPLES)
    TARGET_LINK_LIBRARIES (1elem omega_h)
    TARGET_LINK_LIBRARIES (4elem omega_h)
  endif()
  if(UNIX)
      set(MMG3D_LIB "libmmg3d.so")

  # Windows
  elseif(WIN32)
      set(MMG3D_LIB "mmg3d.dll")
  endif()
    if (REMESH_OMEGA_H)
    TARGET_LINK_LIBRARIES (WeldFormFEM PRIVATE omega_h ${MMG3D_LIB}  OpenMP::OpenMP_CXX)
    else() #ONLY MMG (default)
    TARGET_LINK_LIBRARIES (WeldFormFEM PRIVATE ${MMG3D_LIB} OpenMP::OpenMP_CXX)
    endif()
  endif()

endif()

if (UNIX AND NOT APPLE)
    message(STATUS "Configuring SCOTCH linking for Linux")

    find_library(SCOTCH_LIBRARY scotch REQUIRED)
    find_library(SCOTCHERR_LIBRARY scotcherr REQUIRED)
    find_library(SCOTCHMETIS_LIBRARY scotchmetis REQUIRED)
    find_library(ESMUMPS_LIBRARY esmumps)

    target_compile_definitions(WeldFormFEM PRIVATE USE_SCOTCH)

    target_link_libraries(WeldFormFEM
        PRIVATE
            ${SCOTCH_LIBRARY}
            ${SCOTCHMETIS_LIBRARY}
            ${SCOTCHERR_LIBRARY}
            ${ESMUMPS_LIBRARY} # may be empty if not found, that's ok
    )
endif()


if (BUILD_OMEGA_H)
  if (MSVC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")  # Use the same runtime library as Omega_h
  endif()
  #TARGET_LINK_LIBRARIES (test_omega_h omega_h)
endif()


###################################### COMMIT
set(GIT_COMMIT_HEADER "${CMAKE_BINARY_DIR}/generated/git_commit.h")

# Create directory for generated headers
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated")

# Custom command to generate git_commit.h with the short commit hash
set(GIT_COMMIT_HEADER "${CMAKE_BINARY_DIR}/generated/git_commit.h")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated")

message(STATUS "Version: ${PROJECT_VERSION}")

#add_custom_target(update_git_commit ALL DEPENDS ${GIT_COMMIT_HEADER})
add_custom_target(update_git_commit ALL
  COMMAND ${CMAKE_COMMAND} -DOUT_HEADER=${GIT_COMMIT_HEADER} -DVERSION=${PROJECT_VERSION} -P ${CMAKE_SOURCE_DIR}/cmake/WriteGitHash.cmake
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Forcing regeneration of git_commit.h"
)
add_dependencies(WeldFormFEM update_git_commit)
target_include_directories(WeldFormFEM PRIVATE "${CMAKE_BINARY_DIR}/generated")

#############################################################



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
