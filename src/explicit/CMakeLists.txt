/*************************************************************************/
/*  CMakeLists.txt                                               */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                              */
/*  https://www.opensourcemech.com                                                                */
/*                                                                       */
/*  Copyright (c) 2025-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/


CMAKE_MINIMUM_REQUIRED (VERSION 3.5)

include_directories(${CMAKE_SOURCE_DIR}/lib/LSDynaReader/src)
set (MMG_ROOT "${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/")
link_directories("${MMG_ROOT}/lib")

set(PROJECT_VERSION "0.0.1")  # or your current version
message(STATUS "Version: ${PROJECT_VERSION}")

if (BUILD_GPU)

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

set_source_files_properties(main.C 
                            Solver_explicit.C
                            ${CMAKE_SOURCE_DIR}/include

                            PROPERTIES LANGUAGE CXX)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
ENDIF() 


            
add_executable(weldform_exp main.C
                            Solver_explicit.C
                            ../common/NastranReader.cpp
                         "${SOURCES}"
               )
target_include_directories(weldform_exp PRIVATE "${CMAKE_BINARY_DIR}/generated")

link_directories (${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/lib)


TARGET_LINK_LIBRARIES (weldform_exp PRIVATE 
                                  WeldFormFEM
                                  ${MMG3D_LIB} 
                                  OpenMP::OpenMP_CXX)


add_custom_command(TARGET weldform_exp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:weldform_exp>
        ${CMAKE_BINARY_DIR}/$<TARGET_FILE_NAME:weldform_exp>
    COMMENT "Copying exlicit to root  /"
)
