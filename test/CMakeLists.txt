PROJECT(test)
cmake_minimum_required(VERSION 2.6)
enable_language(CXX Fortran)
# Setting the compilers
# set (CMAKE_Fortran_COMPILER /usr/bin/gfortran)
# set (CMAKE_CXX_COMPILER /usr/bin/g++)
# Setting the flags
set (CMAKE_CXX_FLAGS "-lgfortran")
set_source_files_properties(main.f90 func.cpp PROPERTIES COMPILE_FLAGS -c)
# Making the executable
# ADD_EXECUTABLE(test main.f90 NastranReader.cpp func.cpp addition.cpp addition.h)
# SET_TARGET_PROPERTIES(test PROPERTIES LINKER_LANGUAGE Fortran)

ADD_EXECUTABLE(reader main_reader.cpp NastranReader.cpp )

