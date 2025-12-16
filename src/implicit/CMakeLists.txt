CMAKE_MINIMUM_REQUIRED (VERSION 3.18)

set (BUILD_EXAMPLES OFF)

project (WeldFormImp LANGUAGES CXX )

include_directories(${PROJECT_SOURCE_DIR})

include_directories(../lib/LSDynaReader/src/)

include_directories(../lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h
set (MMG_ROOT "${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/")
link_directories("$MMG_ROOT/lib")
include_directories("${CMAKE_BINARY_DIR}/lib/mmg-5.8.0/include") #FOR: Omega_h_config.h

set(OMEGA_H_DIR "${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src")  # Adjust this path
link_directories(${OMEGA_H_DIR})

set(CMAKE_PREFIX_PATH "/path/to/omega_h")
include_directories(${CMAKE_BINARY_DIR}/lib/omega_h-9.34.13/src) #FOR: Omega_h_config.h

include_directories ("${CMAKE_SOURCE_DIR}/include/common")
include_directories ("${CMAKE_SOURCE_DIR}/src/common")
include_directories ("${CMAKE_SOURCE_DIR}/include/implicit")
include_directories (${CMAKE_SOURCE_DIR}/lib/LSDynaReader/src/)


include_directories(${EIGEN_INC_DIR})
message ("EIGEN DIR" ${EIGEN_INC_DIR})

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
                            Solver.C
                            ${CMAKE_SOURCE_DIR}/include

                            PROPERTIES LANGUAGE CXX)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
ENDIF() 

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
set_source_files_properties(*.C 
                           include/*.h src/include/defs.h src/include/utils.h
                          ${PROJECT_SOURCE_DIR}/src/include/*.cuh

                          PROPERTIES LANGUAGE CXX)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")


option (USE_PETSC "Enable PETSc support (required for parallel)" OFF)
option (USE_METIS "Enable METIS support	                       " OFF)



#This must be included before including another thing else. By now
if (USE_PETSC)
	find_package (PETSc REQUIRED HINTS "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/conf")
	include_directories (..)
	include_directories ("${MY_PETSC_DIR}/include" "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/include" ${PETSC_INCLUDE_PATH})
	
	#THIS IS NOT USED IN WINDOWS, USED MSMPI
	#include_directories ("${MY_PETSC_DIR}/include/mpiuni")
	
	include_directories ("E:/Programas/cygwin64/petsc-3.4.4/include")
	find_library (PETSC_LIB petsc PATH "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib/") 
	if (${PETSC_LIB} STREQUAL "PETSC_LIB-NOTFOUND")
		message (FATAL_ERROR "PETSc library not found")
	endif ()
	##set (EXT_LIBS ${EXT_LIBS} ${PETSC_LIB} ${PETSC_PACKAGE_LIBS})

	#Blas And Lapack, Manual, can be also automatic
	set (PETSC_LIBRARY_DIRS  ${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib)
	link_directories(${PETSC_LIBRARY_DIRS})
	set (EXT_LIBS ${EXT_LIBS} ${PETSC_LIB} flapack fblas)

	set (MODULE_LIST ${MODULE_LIST} "PETSc")
	set (USE_MPI ON)
	add_definitions (-D__PETSC_MODULE)


	link_directories(${MPI_LIB_DIR})
	include_directories(${MPI_INC_DIR})

	set (EXT_LIBS ${EXT_LIBS} msmpi gfortran)
	

endif ()



set (SOURCES                          

                        Matrices.C
                        Solver.C
                        Solver_GlobalMatrix_EPHV.C
                        Solver_GlobalMatrix_QS_martins.C
                        Solver_Eigen.C
                        
                          
                   )
add_library (implicit STATIC "${SOURCES}")                   
add_executable(weldform_imp main.C
                         
               )
               
TARGET_LINK_LIBRARIES (weldform_imp PRIVATE implicit ${MMG3D_LIB} OpenMP::OpenMP_CXX WeldFormFEM)

add_executable(test_1el_imp #ForgeFormFEM.C
                          test_1elem.cpp
               )
TARGET_LINK_LIBRARIES (test_1el_imp PRIVATE implicit ${MMG3D_LIB} OpenMP::OpenMP_CXX WeldFormFEM)

add_executable(test_elastic_cylinder #ForgeFormFEM.C
                          test_elastic_cylinder.cpp
               )
TARGET_LINK_LIBRARIES (test_elastic_cylinder  PRIVATE implicit ${MMG3D_LIB} OpenMP::OpenMP_CXX WeldFormFEM)

add_executable(test_increm_cylinder #ForgeFormFEM.C
                          test_increm_cylinder.cpp
               )
TARGET_LINK_LIBRARIES (test_increm_cylinder  PRIVATE implicit ${MMG3D_LIB} OpenMP::OpenMP_CXX WeldFormFEM)
