/*************************************************************************/
/*  README.md                                                    */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                */
/*  ('https://www.opensourcemech.com',)                                    */
/*                                                                       */
/*  Copyright (c) 2023-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. */ 
/*  See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/




# WelfFormFEM Implicit Version

Updated Lagrangian Implicit Pseudo Static Finite Element Method (FEM) CPU (And Future GPU) based solver. \
WeldFormFEM is aimed to solve solid mechanics large strain problems, such as metal forming. \
The idea is to work via Pure lagrangian with remeshing (adaptive mesh refinement or AMR). \
THIS IS A WORK IN PROGRESS \
WeldFormFEM works both on Linux and Windows.\
I'm building Eigen and PETSC (I have compiled it in both OS).


## Planned Features (WIP)
- Structure Of Arrays (SOA) data arrangement which allows fast CUDA accesing
- Geometric a Material nonlinearities
- Pseudo time to in include materials with strain, strain rate and Temp dependent
- Thermal-Mechanical Coupling
- C++/CUDA CPU/GPU(WIP) Architectures
- Constant Stress Tetra/Triangle Element with Average Nodal Pressure (ANP) for volumetric locking fixing
- Reduced Integration Hexaheadra with viscous hourglass control
- OpenMP (WIP) CPU  parallelization
- Contact algorithm (WIP)
- Automatic Remeshing (rezoning) (WIP)
- 2D Plain Strain/Axisymm & 3D

## HOWTO BUILD 


## CASES WORKING 
Special Cases
test_increm_cylinder.cpp

