# WeldFormFEM
Updated Lagrangian Explicit Finite Element Method (FEM) CPU/GPU based solver. \
WeldFormFEM is aimed to solve solid mechanics large strain problems, such as metal forming. \
The idea is to work via 2 different solvers:
1 - Pure lagrangian with rezoning (adaptive mesh refinement or AMR), Work In Progress as first option \
2 - Coupled Eulerian (fixed mesh or Arbitrarian Eulerian Lagrangian) Lagrangian solver based on Benson works. \
WeldFormFEM works both on Ubuntu and Windows.
you can select to build it to CPU and GPU only by changing a single CMAKE var.

## Features
- Structure Of Arrays (SOA) data arrangement which allows fast CUDA accesing
- Explicit time integration
- C++/CUDA CPU/GPU(WIP) Architectures
- Constant Stress Tetra/Triangle Element with Average Nodal Pressure (ANP) for volumetric locking fixing
- Reduced Integration Hexaheadra with viscous hourglass control
- OpenMP (WIP) CPU  parallelization
- Contact algorithm (WIP)
- Automatic Remeshing (rezoning) (WIP)
- 2D Plain Strain/Axisymm & 3D
- Thermal-Mechanical coupling (WIP)

Locking (left) and fixing (right) tetra \
Compression cylinder \
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/20250117_2.png" width="200" height="200">
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/20250117_1.png" width="200" height="200">

Compression cylinder HEXA w/reduced integration - TETRA \
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/images/hexa_pl_strain.png" width="200" height="200">
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/images/hexa_pressure.png" width="200" height="200">
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/images/tetra_pl_strain.png" width="200" height="200">


Bending \
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/20250117_4.png" width="200" height="100">
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/20250117_3.png" width="200" height="100">

Hexa/Quad Hourglass \
<img src="https://github.com/luchete80/WeldFormFEM/blob/master/20240610_2.png" width="200" height="100">

## Build instructions

CUDACXX=/usr/local/cuda-12.3/bin/nvcc cmake ../WeldFormFEM -DBUILD_GPU=ON

To update libraries (LSDynaReader and Math)

git submodule update --init --recursive



Link here #https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/ to see different architectures. 


Cuurrently working axisymmetric with hourglass for area weight in F90 version

  reduced_int = .True.
  call AddBoxLength(0, V, Lx, Ly, 1.0d0, r, rho, h,reduced_int)

  !! AFTER ADD BOXLEN
  axisymm_vol_weight = .false.
  bind_dom_type = 3 !!!AXISYMM, AFTER CREATING BOX!

  
  elem%sigma_y(:,:) = 300.0e6
  
  do i=1,node_count
  print *, "NODE ELEMENTS "
    print *,"i count ", i , nod%elxnod(i),nod%nodel(i,:)
  end do

  nod%is_fix(:,3) = .true.
 
 
Then in calc elem forces:

'''
if Area weight
#-------------

              fa = 0.25d0/elem%radius(e,gp) * elem%detJ(e,gp) !!! THEN IS WEIGHTED BY 4 in case of gauss point =1
              !!! AREA WEIGHTED, BENSON EQN 2.4.3.2
              !!! 2.4.3.2 remains sig * Area/(4 r0), which is (4detJ)/(4r0) = detJ /r0
              !!! LATER IS MULTIPLIED BY WEIGHT WICH GIVES THE AREA

              elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) - &
                                                     (elem%sigma (e,gp, 1,1) - elem%sigma (e,gp, 3,3) ) * fa
                                                     
              elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) - &
                                                     elem%sigma (e,gp, 1,2) * fa          
              ! print *, "fa ", elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2)
              ! print *, "term 2 ", elem%sigma (e,gp, 1,2) * fa    
              
'''
