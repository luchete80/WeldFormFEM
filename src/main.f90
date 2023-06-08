program WeldFormSPH
!use Integrate
use ElementData
use Domain 
!use Neighbor 
!use Thermal
use omp_lib
use Matrices
use SolverLeapfrog
use VTKOutput
!use SolverVerlet
!use Thermal
!use Mechanical
!use ModPrecision, only : fp_kind

implicit none

  real(fp_kind), dimension(1:3) :: V
  real(fp_kind) :: dx, r, Rxy, Lz, h , ck, young, poisson
  integer:: i, tnr, maxt ,nn
  real(fp_kind),allocatable, dimension(:):: dTdt
  real(fp_kind) :: t_, deltat
  real(fp_kind) :: start, finish
  real(fp_kind) :: L, rho, dt, tf, mat_modK, mat_modG, mat_cs
  
  !type(Domain)	::dom

   ! Variables for clock
   integer count_0, count_1
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
   L = 2.0	
  dx    = 1.0
  r = dx 
  h = dx * 1.2

   V(1) = 0.;V(2) = 0.;V(3) = 0.	
!  !AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)		
  !BOUNDARY CONDITIONS
  !GLOBAL TOP RIGHT NODE , Vx 1m/s, Vy 0.5 m/seconds
  
  rho = 7850.0
  

  
  !PLANE STRESS
  !Plain Stress
  poisson = 0.3
  young = 200.0e9
  allocate (mat_C(3,3))
	ck = young / (1 - poisson*poisson);
	mat_C(1,1) = ck;  mat_C(2,2) = ck;
	mat_C(1,2) = ck*poisson; mat_C(2,1) = ck*poisson;
	mat_C(3,3) = ck*(1.0 - poisson) / 2.0;
  
  print *, "mat_C", mat_C

	! //Plain Strain
	! ck = E*(1. - nu) / ((1. + nu)*(1. - 2.0 * nu));
	! c[0][0] = c[1][1] = ck;
	! c[0][1] = c[1][0] = ck*nu / (1. - nu);
	! c[2][2] = ck*(1. - 2. * nu) / (2.*(1. - nu));
  
  call AddBoxLength(0, V, L, L, L, r, rho, h)
  !!!call AddBoxLength(0, V, L, L, L, r, rho, h)
  
  ! !TODO: CHANGE THIS TO AN ONLY VAULUE, FUNCTION, ETC.
  !CHANGE TO IS_FIXED
  ! nod%is_bcv(1,:) = .true. !Node 1 restricted in 2 dimensions
  ! nod%is_bcv(4,:) = .true. !Node 1 restricted in 2 dimensions

  nod%is_fix(1,:) = .true. !Node 1 restricted in 2 dimensions
  nod%is_fix(2,:) = .true. !Node 1 restricted in 2 dimensions
  nod%is_fix(4,:) = .true. !Node 1 restricted in 2 dimensions
  
  
  ! IF VELOCITY IS APPLIED
  !nod%is_bcv(6,:) = [.false.,.false.,.true.] !GLOBAL DOF TO ADJUST VELOCITY IN A 2 ELEMENT LENGTH CANTILEVDR BEAM  
  !nod%bcv(6,:) = [0.0d0,0.0d0,-1.0d0]
  !fext_glob(6,:) = []
  elem%f_ext(1,6,:) = [0.0d0,0.0d0,-100.0d0]
  
  print *, "BCV 6 ", nod%bcv(6,2)
  print *, "Calculating element matrices "
  

  nod%a(:,:) = 0.0d0
  
  mat_modK= young / ( 3.0*(1.0 -2.0*poisson) );
  mat_modG= young / (2.0* (1.0 + poisson));
  
  mat_cs = sqrt(mat_modK/rho)
  
  elem%cs(:) = mat_cs
  
  dt = 0.7 * dx/(mat_cs)
  tf = dt * 3.0
  
  elem%rho(:) = rho
  
  print *, "Shear and Bulk modulus", mat_modG,mat_modK
  print *, "time step size with CFL 0.7", dt
  call SolveLeapfrog(tf,dt)
  
  call WriteMeshVTU
  
  open (1,file='test.csv')!, position='APPEND')  
  write (1,*) "X, Y, Z"
  
  do i=1,node_count
    write (1,*) nod%x(i,1), ", ", nod%x(i,2), ", " ,nod%x(i,3)  
    end do
  close(1)

!  do i = 1, part_count
!  !write (*,*) "Particle", i ," position is ", pt%x(i,1), pt%x(i,1), pt%x(i,3)
!  end do 
!  !call AddCylinderLength(0, V, Rxy, Lz, r)

!  call DomInit(12)
!  call InitNb()
!  call AllocateNbData()
!  
!  call CellInitiate()
!  call ListGenerate()
!  
!  call MainNeighbourSearch()
!  call InitRedArraysOnce()
!  call CalcPairPosList()
!  
!  allocate (dTdt(part_count))
!  allocate (temp_pair(pair_tot_count))
!  pt%t(:)     = 20.
!  pt%cp_t(:)  = 1.
!  pt%k_t(:)   = 3000.  
!  
!  print *, "Size of floating point: ", sizeof(pt%t(1))

!  !call cpu_time(start)
!  !deltat = 0.00036
!  deltat = 0.3*h*h*rho*pt%cp_t(1)/pt%k_t(1)	
!  
!  print *,'Reduction by particle... '    
!  t_ = 0.

!   ! Starting time
!   call system_clock(count_0, count_rate, count_max)
!   time_init = count_0*1.0/count_rate

!  
!  nn=(L/dx)*(L/dx)
!    pt%t(1:nn) = 500.
!  !!! FASTEST ALGORITHM
!  do while (t_ <= 1000.0*deltat)

!    call CalcTempIncPart(dTdt) !THIS IS THE FASTEST BY NOW
!    !print *, "dTdt 0",  dTdt(1)
!    pt%t(:) = pt%t(:) + dTdt(:) * deltat
!    
!    t_ = t_ + deltat
!  end do

!  ! Ending time
!  call system_clock(count_1, count_rate, count_max)
!  time_final = count_1*1.0/count_rate
!  ! Elapsed time
!  elapsed_time = time_final - time_init

!  ! Write elapsed time
!  write(*,1003) int(elapsed_time),elapsed_time-int(elapsed_time)
!  !write(*,*) "Elapsed time: ", int(elapsed_time),elapsed_time-int(elapsed_time)
!   
!  !call cpu_time(finish)
!  print *,'Time: ', t_ 
!  print '("CPU Time = ",f6.3," seconds.")',finish-start
!  print *, "Program End."

!  
!  open (1,file='temp.csv')!, position='APPEND')  
!  write (1,*) "X, Y, Z, temp"

!  do i=1,part_count  
!    write (1,*) pt%x(i,1), ", ", pt%x(i,2), ", " ,pt%x(i,3), ", " ,pt%t(i) 
!  end do
!  close(1) 



!  1003 format('  Elapsed Time  = ',i0,f0.9)
!  
end program WeldFormSPH
