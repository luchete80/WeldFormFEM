!!!!BENSON 1992

! (1) Knowing the stress, pressure, hourglass forces and shock viscosity at t” in each zone or
! element, the forces at the nodes are calculated. The accelerations of the nodes are
! calculated by dividing the nodal forces by the nodal masses.
! (2) The acceleration is integrated to give the velocity at t”+l/2”.
! (3) The velocity is integrated to give the displacement at t”+‘.
! (4) The constitutive model for the strength of the material is integrated from t to t_n+1 now
! that the motion of the material is known.
! (5) The artificial shock viscosity and hourglass viscosity are calculated from un+1/2. ATTENTION
! (6) The internal energy is updated based on the work done between tn and t_n+1.
! (7) Based on the density and energy at t_n+l, the pressure is calculated from the equation of
! state.
! (8) A new time step size is calculated based on the speed of sound through each of the
! elements and their geometry.
! (9) Advance the time and return to step (1)

!! VELOCITY VERLET
!!Calculate v(t+1/2dt) = v(t) +1/2a(t) dt
!!Calc      x(t+dt ) = x(t) + v(t+1/2dt) dt
!!obtain    a(t+dt) using x(t+dt)
!!update    v(t+dt) = v(t+1/2dt)+1/2a(t+dt) dt


module SolverVerlet
use ModPrecision, only : fp_kind

contains 

subroutine SolveVerlet (tf, dt)
  use omp_lib
  use Matrices
  use Mechanical
  
  implicit none
  integer :: n, d, iglob
  
  real(fp_kind),intent(in)::tf, dt
  
  real(fp_kind), dimension(node_count*dim) :: mdiag !!DIAGONALIZATION COULD BE DONE INSIDE ACC CALC  
  real(fp_kind), dimension(dim) :: prev_acc
 

  call set_edof_from_elnod()
  
  call calculate_element_Jacobian()
  call calculate_element_shapeMat() !AND MASS
  
  ! call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU
  call assemble_mass_matrix()
    mdiag(:)=0.0d0
    do iglob =1, (node_count * dim)
      n = 1
      do while (n .le. node_count * dim) !column
         mdiag(iglob) = mdiag(iglob) + m_glob(iglob,n)
         n = n+ 1
      end do !col
    end do   
  calc_m = .False.
 !print *, "M Diag ", mdiag
  !print *, "m glob", m_glob
  
  nod%u(:,:) = 0.0d0
  
  time = 0.0  
  print *,"main loop"
  do while (time .le. tf)

    call calculate_element_Jacobian()
    call calculate_element_derivMat()
    ! call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU
    ! !NODAL CALCULATION

    ! call assemble_int_forces()
  ! (1) Knowing the stress, pressure, hourglass forces and shock viscosity at t” in each zone or
  ! element, the forces at the nodes are calculated. The accelerations of the nodes are
  ! calculated by dividing the nodal forces by the nodal masses.   
    print *, "calc elem forces "
    call cal_elem_forces()
    print *, "assemble int forces "
    call assemble_int_forces()
    
    print *, "calc accel "
    do n=1,node_count
      !prev_acc(:) = nod%a(n,:)
      do d=1,dim
        iglob = (n-1) * dim + d !TODO: CALL THIS IN A FUNCTION
        nod%a(n,d) = rint_glob(n,d)/mdiag(iglob) 
      end do 
    end do

  !!!!! THIS IS NOT SOLVED AS A COMPLETED STEP (REDUCED VERLET=
  ! (2) The acceleration is integrated to give the velocity at t”+l/2”.
    ! !Update vel with CURRENT ACCELERATION
    ! THIS WOULD BE AT ONE STEP
    ! nod%v(n,:) = nod%v(n,:) + dt * 0.5 * (nod%a(n,:) + prev_acc(:)) 
    ! print *,"node vel ", nod%v(n,:)  
    nod%v(n,:) = nod%v(n,:) + dt * 0.5 * nod%a(n,:)

  !!(3) The velocity is integrated to give the displacement at tn+1.
  nod%x(n,:) = nod%x(n,:) +  nod%v(n,:) * dt
  
  call calc_elem_vol()
  call disassemble_uvele()     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains ()     !!!!!STRAIN AND STRAIN RATES
  
  ! (4) The constitutive model for the strength of the material is integrated from t to t_n+1 now
  ! that the motion of the material is known.
  
  ! (5) The artificial shock viscosity and hourglass viscosity are calculated from un+1/2. ATTENTION
  call calc_hourglass_forces
  
    ! call assemble_int_forces()
    ! ! !Diagonalize
    ! ! !SIMPLEST FORM, ROW SUM 
    ! !print *, "Calc mdiag"
    ! ! if (calc_m .eqv. .True.) then 
      ! ! call assemble_mass_matrix()
      ! ! mdiag(:)=0.0d0
      ! ! iglob = 1
      ! ! do while (iglob .le. node_count * dim)
        ! ! n = 1
        ! ! do while (n .le. node_count * dim) !column
           ! ! mdiag(iglob) = mdiag(iglob) + m_glob(iglob,n)
           ! ! n = n+ 1
        ! ! end do !col
      ! ! iglob = iglob + 1
      ! ! end do 
    ! ! end if
    
    ! print *, " mglob ," ,m_glob
    ! print *, " mdiag ," ,mdiag
    
    ! call impose_bcv
    ! !Calculate positions with PREVIOUS ACCELERATIONS
    ! n = 1
    ! print *, "calculating positions "
    ! do while (n .le. node_count)
      ! !print *, "n ", n, "acc ", nod%v(n,:)
      ! nod%u(n,:) = nod%u(n,:) + nod%v(n,:) * dt + nod%a(n,:) * dt * dt  
      ! nod%x(n,:) = nod%x(n,:) + nod%u(n,:)    
      ! !nod%x(n,:) = nod%x(n,:) +  nod%v(n,:) * dt + nod%a(n,:) * dt * dt    
      ! print *,"node u ", nod%u(n,:)
      ! n = n + 1
    ! end do !Node 
    
    ! !Again to calculate strains
    ! call disassemble_uele()
    ! print *, "u ele 2", elem%uele(2,:,:)
    ! call cal_elem_strains ()
    
    ! print *,"strains", elem%strain(:,:,:,:)
    ! !Calculate Nodal accelerations a(t+dt) from rext(t)-rint(t)-fcont
    ! print *, "calculating accel"
    
    ! do n=1,node_count
      ! !print *, "node ", n 
      ! prev_acc(:) = nod%a(n,:)
      ! d = 1
      ! do while (d .le. 2)
        ! !print *, "dim ", d
        ! iglob = (n-1) * dim + d !TODO: CALL THIS IN A FUNCTION
        ! print *, "iglob ", iglob, "mdiag ", mdiag(iglob)
        
        ! nod%a(n,d) = rint_glob(n,d)/mdiag(iglob) 
        ! print *,"node acc ", nod%a(n,:), "rint ", rint_glob(n,:)
        ! d = d + 1
      ! end do
      ! !Solve motion at t_n+1
      ! !Update vel with CURRENT ACCELERATION
      ! nod%v(n,:) = nod%v(n,:) + dt * 0.5 * (nod%a(n,:) + prev_acc(:)) 
      ! print *,"node vel ", nod%v(n,:)

    ! end do !Node
    ! call impose_bca
    
    ! !REINFORCE bc velocity
    ! !TODO: SEPARATE FUNCTION
    ! call impose_bcv
    
    ! !calc_m = .False.
  time = time + dt
  end do !time

end subroutine SolveVerlet

end module SolverVerlet