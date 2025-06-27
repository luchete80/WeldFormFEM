# ForgeFormFEM

Updated Lagrangian Implicit Pseudo Static Finite Element Method (FEM) CPU (And Future GPU) based solver. \
WeldFormFEM is aimed to solve solid mechanics large strain problems, such as metal forming. \
The idea is to work via Pure lagrangian with remeshing (adaptive mesh refinement or AMR). \
THIS IS A WORK IN PROGRESS \
WeldFormFEM works both on Ubuntu and Windows.\
I'm building Eigen and PETSC (I have compiled it in both OS).


## Planned Features (WIP)
- Structure Of Arrays (SOA) data arrangement which allows fast CUDA accesing
- Geometric a Material nonlinearities
- Pseudo time to in include materials with strain, strain rate and Temp dependent
- Thermal Coupled 
- C++/CUDA CPU/GPU(WIP) Architectures
- Constant Stress Tetra/Triangle Element with Average Nodal Pressure (ANP) for volumetric locking fixing
- Reduced Integration Hexaheadra with viscous hourglass control
- OpenMP (WIP) CPU  parallelization
- Contact algorithm (WIP)
- Automatic Remeshing (rezoning) (WIP)
- 2D Plain Strain/Axisymm & 3D
- Thermal-Mechanical coupling (WIP)

## HOWTO BUILD 

