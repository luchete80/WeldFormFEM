project (examples LANGUAGES CXX )
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
include_directories(${PROJECT_SOURCE_DIR}/src/include)
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
