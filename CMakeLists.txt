cmake_minimum_required(VERSION 3.2)

project(WeldFormSPH)
enable_language(Fortran CXX)
                            
                            
INCLUDE      (FindOpenMP)

if(OPENMP_FOUND)
    ADD_DEFINITIONS (-DUSE_OMP)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
        SET (MISSING "${MISSING} OpenMP")
endif(OPENMP_FOUND)

#if GNU, 
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall"})
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
endif()

file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(WeldForm ${sources})

SET_TARGET_PROPERTIES(WeldForm PROPERTIES LINKER_LANGUAGE Fortran)