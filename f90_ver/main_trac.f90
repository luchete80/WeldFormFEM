program WeldFormFEM
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
use class_ContMesh
use Mechanical
use SolverChungHulbert
!use SolverVerlet
!use Thermal
!use Mechanical
!use ModPrecision, only : fp_kind

implicit none

   interface

      SUBROUTINE reader(fname, nod, nodecount, elnod) BIND(C, name="ReadNastranTriMesh") ! ****
      import :: c_ptr, c_char
      character(C_CHAR), intent(in)  :: fname(*)
      TYPE(c_ptr) , intent(out):: nod ! ****
      type(c_ptr), intent(out)::nodecount
      TYPE(c_ptr) , intent(out):: elnod
      !logical, intent(in) :: issurf
      !type(c_int), intent(in) :: dimm
      

      END SUBROUTINE reader
      
   end interface
  
  real(fp_kind), pointer :: nodetest(:)
  integer, pointer :: eltest(:)
  real(fp_kind), dimension(1:3) :: V
  real(fp_kind) :: dx, r, Rxy, h , ck, young, poisson
  integer:: i, tnr, maxt ,nn
  real(fp_kind),allocatable, dimension(:):: dTdt
  real(fp_kind) :: t_, deltat
  real(fp_kind) :: start, finish
  real(fp_kind) :: Lx, Ly, Lz, rho, dt, tf, mat_modG, mat_cs
  logical :: reduced_int
  type(c_ptr) :: nodptr,elnodptr,ncount ! ****
  integer, pointer :: ncount_int
  type(dom_type)	::dom

   ! Variables for clock
   integer count_0, count_1, gp
   integer count_rate, count_max
   double precision time_init, time_final, elapsed_time

    
  call omp_set_num_threads(12);
  
  maxt = omp_get_max_threads()
  
   write( *, * ) 'Max threads ', maxt
   
  !CALL 10 x 10 element rectangle
  ! dx    = 0.1
  ! Rxy  = 0.15
  ! Lz = 0.56
  ! r = dx / 2.0
  ! h = dx * 1.2
  !!!! 2 ELEMENT LENGTH CANTILEVDR BEAM

  Dim = 3
  Lx = 0.01	
  Ly = 0.01
  Lz = 0.1
  dx    = 0.012d0
  r = dx /2.0
  h = dx * 1.2

   V(1) = 0.;V(2) = 0.;V(3) = 0.	
!  !AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)		
  !BOUNDARY CONDITIONS
  !GLOBAL TOP RIGHT NODE , Vx 1m/s, Vy 0.5 m/seconds
  
  rho = 2700.0
  
  poisson = 0.3
  young = 70.0e9
  
  print *, "mat_C", mat_C

	! //Plain Strain
	! ck = E*(1. - nu) / ((1. + nu)*(1. - 2.0 * nu));
	! c[0][0] = c[1][1] = ck;
	! c[0][1] = c[1][0] = ck*nu / (1. - nu);
	! c[2][2] = ck*(1. - 2. * nu) / (2.*(1. - nu));
  
  reduced_int = .false.
  call AddBoxLength(0, V, Lx, Ly, Lz, r, rho, h,reduced_int)
  
  do i=1,node_count
  print *, "NODE ELEMENTS "
    print *,"i count ", i , nod%elxnod(i),nod%nodel(i,:)
  end do

  do i=1,node_count
    if (nod%x(i,3)<r) then
      nod%is_fix(i,:) = .true.
      print *, "Fixed Node ", i, " at: ", nod%x(i,:)
    end if
  end do

  do i=1,node_count
    if (nod%x(i,3)> (Lz - r) ) then
      nod%is_bcv(i,3) = .true.
      nod%bcv(i,3) = 0.1d0
      print *, "Velocity Node ", i, " at: ", nod%x(i,:)
    end if
  end do
  
  ! nod%bcv(3,:) = [0.0d0,-1.0d0]
  ! nod%bcv(4,:) = [0.0d0,-1.0d0]
  
  ! nod%is_fix(1,:) = .true. !Node 1 restricted in 2 dimensions
  ! nod%is_fix(2,2) = .true. !Node 1 restricted in 2 dimensions
  
  
 ! print *, "BCV 6 ", nod%bcv(6,3)
  print *, "Calculating element matrices "
  

  nod%a(:,:) = 0.0d0
  
  mat_K= young / ( 3.0*(1.0 -2.0*poisson) );
  mat_G= young / (2.0* (1.0 + poisson));
  
  mat_cs = sqrt(mat_K/rho)
  mat_cs0 = mat_cs
  print *, "Material Cs: ", mat_cs
  
  elem%cs(:) = mat_cs
  
  !dt = 0.7 * dx/(mat_cs)
  dt = 0.7 * dx/(mat_cs)
  tf = dt * 1.0
  
  tf = 1.0e-7
  tf = 1.0e-5
  
  elem%rho(:,:) = rho
  
  print *, "Shear and Bulk modulus", mat_modG,mat_K
  print *, "time step size with CFL 0.7", dt
  !call SolveLeapfrog(tf,dt)
  !call SolveVerlet(dom,tf,dt)
  call SolveChungHulbert(tf,dt)
  call CalcEquivalentStress()
  call AverageData(elem%rho(:,1),nod%rho(:))
  call AverageData(elem%sigma_eq(:,1),nod%sigma_eq(:))
    
  call WriteMeshVTU('output.vtu')
  
  open (1,file='test.csv')!, position='APPEND')  
  write (1,*) "X, Y, Z"
  
  do i=1,node_count
    write (1,*) nod%x(i,1), ", ", nod%x(i,2), ", " ,nod%x(i,3)  
    end do
  close(1)
  
  print *, "Element stresses"
  do i=1,elem_count
    do gp=1, elem%gausspc(i)
      print *, elem%sigma(i,gp,:,:)
    end do
  end do

  do i=1,node_count
    print *, "nod ", i, "Disp ", nod%u(i,:)  
  end do  
  
 
end program WeldFormFEM
