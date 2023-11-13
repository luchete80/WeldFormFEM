program contact
  use iso_c_binding
  
!use Integrate
use ElementData
use Domain 
!use Neighbor 
!use Thermal
use omp_lib
use Matrices
use SolverLeapfrog
use SolverVerlet
use VTKOutput
!use class_ContMesh
use Mechanical  !Calc equivalent stress
use class_ContMesh
!use SolverVerlet
!use Thermal
!use Mechanical
!use ModPrecision, only : fp_kind

implicit none

  type(Mesh) :: msh
  integer::axis=3
  logical:: positaxisorent = .true.
  real(fp_kind),dimension(3) :: p1 = [0,0,0], p2 = [1,1,0]
  integer :: dens = 2
  
  Dim =2 
  
  call  AxisPlaneMesh(msh, axis, positaxisorent, p1, p2,  dens)
	
	print *, "Contact elements ", msh%elem_count, ", nodes: ", msh%node_count
 
  print *, "End program. "
  
end program contact