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
use SolverChungHulbert
use SolverGreenNag
use SolverKickDrift
!use SolverRedVerlet
use VTKOutput
!use class_ContMesh
use Mechanical  !Calc equivalent stress
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
  real(fp_kind) :: dx, r, Rxy, Lz, h , ck, young, poisson
  integer:: i, tnr, maxt ,nn, el, d, d2
  real(fp_kind),allocatable, dimension(:):: dTdt
  real(fp_kind) :: t_, deltat
  real(fp_kind) :: start, finish, sy
  real(fp_kind) :: L, Rtot, rho, dt, tf, mat_modK, mat_modG, mat_cs
  logical :: reduced_int, dim_2D
  type(c_ptr) :: nodptr,elnodptr,ncount ! ****
  integer, pointer :: ncount_int  
  integer :: nc_1, nc_2, nc_3
  !type(Domain)	::dom

   ! Variables for clock
   integer count_0, count_1, gp
   integer count_rate, count_max
   integer nc1, nc2, nc3
   double precision time_init, time_final, elapsed_time

  type(Dom_type) :: dom
  
  call omp_set_num_threads(12);
  
  maxt = omp_get_max_threads()
	
  dim_2D = .True.
   write( *, * ) 'Max threads ', maxt
   
  !CALL 10 x 10 element rectangle
  ! dx    = 0.1
  ! Rxy  = 0.15
  ! Lz = 0.56
  ! r = dx / 2.0
  ! h = dx * 1.2
  !!!! 2 ELEMENT LENGTH CANTILEVDR BEAM

  Dim = 2
  L = 0.616
  Rtot = 0.15
  
  !dx    = 0.05d0
  dx    = 0.075d0
  r = dx /2.0

  V(1) = 0.;V(2) = 0.;V(3) = 0.	
!  !AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)		
  !BOUNDARY CONDITIONS
  !GLOBAL TOP RIGHT NODE , Vx 1m/s, Vy 0.5 m/seconds
  
  rho = 2700.0
  
  poisson = 0.3
  young = 206.0e9
  sy = 300.0e6
  
  print *, "mat_C", mat_C


	! //Plain Strain
	! ck = E*(1. - nu) / ((1. + nu)*(1. - 2.0 * nu));
	! c[0][0] = c[1][1] = ck;
	! c[0][1] = c[1][0] = ck*nu / (1. - nu);
	! c[2][2] = ck*(1. - 2. * nu) / (2.*(1. - nu));
  
  reduced_int = .true.
  !!!! DO NOT PUT z = 0.0d0, it crashes
  call AddBoxLength(0, V, Rtot, L, 1.0d0, r, rho, h,reduced_int)
  bind_dom_type = 3 !!!AXISYMM, AFTER DOMAIN CREATION!!!
    
  nc1 = 0;nc2 = 0;nc3 = 0;
  print *, "NODE ELEMENTS "
    do i=1,node_count
    ! print *,"i count ", i , nod%elxnod(i),nod%nodel(i,:)
    
    if (nod%x(i,1) < r) then     
      nod%is_fix(i,1) = .true. !AXI SYMMETRIC
      nc1 = nc1 + 1
    end if
        
    if (nod%x(i,2) < r) then 
      nod%is_fix(i,:) = .true. !Node 1 restricted in 2 dimensions, AXIS AND VERTICAL  
      nc2 = nc2 +1      
    end if

    if (nod%x(i,2) > L - r) then 
      nod%bcv(i,:) = [0.0d0,-1.0d0]
      nod%is_bcv(i,2) = .true.
      nod%is_fix(i,1) = .true. 
      nc3 = nc3 + 1 
    end if
    
  end do
  
  print *, "Boundary nodes ", nc1, ", ", nc2, ", ", nc3
    

  
 ! print *, "BCV 6 ", nod%bcv(6,3)
  print *, "Calculating element matrices "
  

  nod%a(:,:) = 0.0d0
  
  !ONLY FOR ELASTIC
  dom%mat_E = young
  dom%mat_nu = poisson
  
  mat_modK= young / ( 3.0*(1.0 -2.0*poisson) );
  mat_G= young / (2.0* (1.0 + poisson));
  
  dom%mat_K = mat_modK !!!TODO CREATE MATERIAL
  
  mat_cs = sqrt(mat_modK/rho)
  mat_cs0 = mat_cs
  print *, "Material Cs: ", mat_cs
  
  elem%cs(:) = mat_cs
  
  dt = 0.3 * dx/(mat_cs)
  print *, "time step size with CFL 0.7", dt
  
  !dt = 5.0e-6
  !tf = 1.5e-4

  !dt = 1.0e-5
  !tf = 1.0e-3 
  tf = 2.0*dt 
  
  elem%rho(:,:) = rho
  
  print *, "Shear and Bulk modulus", mat_modG,mat_modK


  !call SolveVerlet(dom,tf,dt)
  !call SolveKickDrift(tf,dt) !!!!!OLD; NOT WORKING
  !call SolveLeapfrog(dom,tf,dt)

  call SolveChungHulbert(dom,tf,dt)
  !call SolveLeapFrog(dom,tf,dt)
  !call SolveGreenNag(dom,tf,dt)
  call CalcEquivalentStress()
  call AverageData(elem%rho(:,1),nod%rho(:))

  do d=1,dim
    call AssembleElNodData(elem%hourg_nodf(:,:,d),nod%f_hour(:,d))
  end do
  
  call AverageData(elem%sigma_eq(:,1),nod%sigma_eq(:))
  
    do d=1,3
      do d2=1,3
      call AverageData(elem%sigma(:,1,d,d2),nod%sigma(:,d,d2))
      end do
    end do
  
  call WriteMeshVTU('output.vtu')
  
  open (1,file='test.csv')!, position='APPEND')  
  write (1,*) "X, Y, Z"

  do i=1,node_count
    print *, "nod ", i, "Disp ", nod%u(i,:)  
  end do  

  do i=1,node_count
    print *, "nod ", i, "Vel ", nod%v(i,:)  
  end do  

  do i=1,node_count
    print *, "nod ", i, "Acc ", nod%a(i,:)  
  end do  
  
  do i=1,node_count
    write (1,*) nod%x(i,1), ", ", nod%x(i,2), ", " ,nod%x(i,3)  
    end do
  close(1)
  
  print *, "Element stresses"
  do i=1,elem_count
    do gp=1, elem%gausspc(i)
      print *, elem%sigma(i,gp,:,:)
      !print *, "Sigma eq ", elem%sigma_eq(i,gp)
    end do
  end do
  print *, "Element strain rates" 
  do i=1,elem_count
    do gp=1, elem%gausspc(i)
      print *, elem%str_rate(i,gp,:,:)
    end do
  end do
  print *, "Accels"
  do i=1,node_count
    print *, "a ", nod%v(i,:)  
  end do  
  
  ! print *, "Element rot rates" 
  ! do i=1,elem_count
    ! do gp=1, elem%gausspc(i)
      ! print *, elem%rot_rate(i,gp,:,:)
    ! end do
  ! end do
  
  print *, "Global forces "
    do nn=1,node_count
        print *, rint_glob(nn,:)
    end do
        
  print *, "Internal forces " 
  do i=1,elem_count  
    print *, "elem ", i
    do nn=1,nodxelem
      print *, elem%f_int(i,nn,:) 
    end do
  end do
  
  ! print *, "Hourglass forces " 
  ! do i=1,elem_count  
    ! print *, "elem ", i
    ! do nn=1,nodxelem
      ! print *, elem%hourg_nodf(i,nn,:) 
    ! end do
  ! end do
  !(fname, node, elnod, dimm, issurf)
  !print *, "dim: ", dim, "is surf "
  allocate (ncount_int)
  !allocate (ncount)
  ! call reader('cylinder.nas', nodptr, elnodptr, ncount)
  ! !print *, "node count ", ncount
  ! CALL C_F_POINTER(ncount, ncount_int)
  
  ! !print *, "Size of ptr", size(nodptr)
  ! CALL C_F_POINTER(nodptr, nodetest, [100])
  ! CALL C_F_POINTER(elnodptr, eltest, [100])
  ! print *, "node count ", ncount_int
  
  !call MeshCSVreader()
  !print *, "Nodetest " , nodetest(1), ", ", nodetest(2), ", ", nodetest(3)
  
!  do i = 1, part_count
!  !write (*,*) "Particle", i ," position is ", pt%x(i,1), pt%x(i,1), pt%x(i,3)
!  end do 
!  !call AddCylinderLength(0, V, Rxy, Lz, r)
 
end program WeldFormFEM
