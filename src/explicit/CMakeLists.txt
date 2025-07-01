CMAKE_MINIMUM_REQUIRED (VERSION 3.5)


include_directories(${CMAKE_SOURCE_DIR}/lib/LSDynaReader/src)

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

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
ENDIF() 


            
add_executable(weldform_exp main.C
                            Solver_explicit
                         "${SOURCES}"
               )
target_include_directories(weldform_exp PRIVATE "${CMAKE_BINARY_DIR}/generated")

TARGET_LINK_LIBRARIES (weldform_exp PRIVATE ${MMG3D_LIB} OpenMP::OpenMP_CXX WeldFormFEM)



