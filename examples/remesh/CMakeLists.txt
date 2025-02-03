project (remesh_ex1 LANGUAGES CXX )
include_directories(include_directories ../../)
include_directories(include_directories ../../lib)
include_directories(include_directories ../../lib/omega_h-9.34.13/src)
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(../../src/include)
include_directories(../../src)

if (BUILD_GPU)
set(CMAKE_CXX_STANDARD 11) #LIKE IN CUNSearch, this allows to include string (conflicts with CXX17)!!!
set(CMAKE_CXX_STANDARD_REQUIRED ON)
else()
set_source_files_properties(*.C ../../src/*.C PROPERTIES LANGUAGE CXX) #SINCE ARE DEFINED AS .C
endif()

add_executable(remesh_ex1
                       remesh_ex.C
                       ../../src/Domain_d.C
                       ../../src/ReMesher.C
                       ../../src/VTKWriter.C
             )

#TARGET_LINK_LIBRARIES (remesh_ex1  omega_h )

if (BUILD_GPU)
#THIS SHOULD BE SET GREATER THAN 60 DUE TO CUDA ADDATOMIC
#set_property(TARGET 4elem PROPERTY CUDA_ARCHITECTURES 75)


#TARGET_LINK_LIBRARIES (remesh_ex1  omega_h cuda cudadevrt)

# TARGET_LINK_LIBRARIES (WeldFormFEM cuda cudadevrt)
else()
# TARGET_LINK_LIBRARIES (MetFEMGPU libstdc++)
endif()
