  ! ! // Compute the Jacobian
  ! ! model->computeJacobian(true);
  ! ! model->computeUnderJacobian(true);
  ! ! model->computeMassMatrix();
  ! ! computeTimeStep(true);
  ! ! model->computeStrains();
  ! ! model->computePressure();
  ! ! model->computeStress(timeStep);
  ! ! model->computeFinalRotation();
  ! ! model->computeInternalForces();
!!! FROM DYNAELA
!!!! :____ MAIN LOOP 
    ! computePredictions();
    ! model->computeStrains(); !With deformation gradient !F.polarCuppenLnU(_integrationPoint->StrainInc, _integrationPoint->R);
    ! model->computePressure(); !Stress.thirdTrace() + K * pressureIncrement; pressureIncrement += StrainInc.trace();
    ! model->computeStress(timeStep);
    ! model->computeFinalRotation();
    ! model->computeInternalForces();
    ! explicitSolve();
    ! computeDensity();


    ! if (model->currentTime < _solveUpToTime)
    ! {

      ! model->computeJacobian();
      ! model->computeUnderJacobian();

      ! computeTimeStep();

    ! }

module SolverGreenNag
use ModPrecision, only : fp_kind

contains 

!!!COMPUTE STEP SIZE
!!!maximumFrequency = 2.0 / model->computeCourantTimeStep();
!!!timeStep = _timeStepSafetyFactor * _omegaS / maximumFrequency;


subroutine SolveGreenNag (domi, tf, dt)
  use omp_lib
  use Matrices
  use Mechanical
  
  implicit none
  integer :: n, d, iglob, step, e, gp
  
  logical :: first_step, x_at_midtime
  logical :: debug_mode 
  real(fp_kind),intent(in)::tf, dt
  
  real(fp_kind), dimension(node_count) :: mdiag !!DIAGONALIZATION COULD BE DONE INSIDE ACC CALC  
  real(fp_kind), dimension(dim) :: prev_acc
  real(fp_kind), dimension(node_count,dim) :: u, prev_a, x_temp
  
  real(fp_kind) :: alpha, beta, gamma, rho_b, omega !!! CHUNG HULBERT PARAMETERS
 
  real(fp_kind), dimension(nodxelem,dim) :: xtest

  type (dom_type), intent (in) :: domi

  call set_edof_from_elnod()
  
  call calculate_element_Jacobian
  print *, "shape mat"
  !call calculate_element_shapeMat() !AND MASS
  call calc_elem_vol !!!! In order to define initial volume
  call calculate_element_derivMat()
  elem%vol_0(:) = elem%vol(:)
  !print *,"Element Initial Vol"
  ! do n = 1, elem_count
    ! print *, elem%vol(n)
  ! end do    
    
  ! call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU
  print *, "Assemblying mass matrix" 
  call assemble_mass_matrix()
  !print * , "done"
  !print *, "mass matrix",m_glob
    mdiag(:)=0.0d0
    do iglob =1, node_count
      do n=1, node_count  !column
         mdiag(iglob) = mdiag(iglob) + m_glob(iglob,n)
      end do !col
    end do   
  calc_m = .False.
  
  !!!! ONLY FOR TESTING
  do n=1, node_count  !column
     mdiag(n) = tot_mass/node_count 
  end do
  
  print *, "M Diag ", mdiag
  
  !print *, "m glob", m_glob
  ! print *, "done"
  nod%u(:,:) = 0.0d0
  debug_mode = .false.
  first_step  = .true.
  
  !!!!!!!!!!!!!!! IF EXTERNAL FORCES (AND IF NOT?????, IF BCs ARE ONLY VELOCITY??
  !!!!!!!!!!!!!! CALCULATE Ku0 = RINT0, Initial internal forces
  ! print *, "Assemblying forces..."
  ! call assemble_forces()
  ! do n=1,node_count
      ! nod%a(n,:) = (fext_glob(n,:)-rint_glob(n,:))/mdiag(n) 
      ! print *, "fext n ", n, fext_glob(n,:)
  ! end do
  ! call impose_bca
  
  ! do n=1,node_count
    ! print *, "Initial accel ", n, "a ", nod%a(n,:)  
  ! end do  
  
  !!!! IF ONLY ARE SET bcv 
  call impose_bcv

  prev_a = nod%a
  
  !!!! IS THERE ANY STRESS?
  elem%sigma (:,:,:,:) = 0.0d0 !!!! FOR INT FORCES (elem%f_int(e,gp,d,d)) CALCULATION
  !elem%f_int (:,:,:)   = 0.0d0 !!!! I Ncal_elem_forces

  rho_b = 0.8182  !!! DEFAULT SPECTRAL RADIUS
  
  alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b)
  beta = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b))
  gamma = 1.5 - alpha;
  omega = sqrt((12.0 * (1.0 + rho_b)**3.0 * (2.0 - rho_b)) / &
                (10.0 + 15.0 * rho_b - rho_b * rho_b + rho_b**3.0 - rho_b**4.0))

  print *, "alpha ", alpha
  print *, "beta ", beta
  print *, "gamma ", gamma
  print *, "omega ", omega
  
  !In central difference (leapfrog)
  
  
  elem%shear_stress = 0.0d0 
  time = 0.0  
  step = 0
  
  x_at_midtime = .False.
  
  print *,"------------------------------------------------------------------------------------------------"
  print *,"main loop,  -------------------------------------------------------------------------------"
  do while (time < tf)
    step = step + 1
    print *, "Time: ", time, ", step: ",step, "---------------------------------------------------------"
    ! if (step == 1) then
      ! dt = tf - time
    ! end if
  ! if (time < 100.0d0*dt) then
    ! nod%bcv(5:8,3) = -1.0 * time/(10.0d0*dt)
    ! !nod%bcv(3:4,2) = -0.1 * time/(100.0d0*dt)
  ! else 
    ! nod%bcv(5:8,3) = -1.0
    ! !nod%bcv(3:4,2) = -0.1
  ! end if 
  
    do n=1,elem_count
      if (elem%gausspc(n) .eq. 8) then !!!! ELSE IS CONSTANT
        call calculate_element_shapeMat() !AND MASS
      end if
    end do
  !!! PREDICTION PHASE
  u = dt * (nod%v + (0.5d0 - beta) * dt * prev_a)
  nod%u = nod%u + u

  print *, "Initial a"
  do n=1,node_count
  print *, prev_a(n,:)
  end do
  
  do n=1,node_count
  print *, "Initial u ", u(n,:)
  end do
  !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
  print *, "Initial v"
  do n=1,node_count
  print *, nod%v(n,:)
  end do
  nod%v = nod%v + (1.0d0-gamma)* dt * prev_a
  nod%a = 0.0d0
  
  call impose_bcv !!!REINFORCE VELOCITY BC
  print *, "Pred v"
  do n=1,node_count
  print *, nod%v(n,:)
  end do
  ! nod%u = nod%u +  nod%v * dt!/2.0  
  ! nod%x = nod%x + nod%u             !! EVALUATE dHdxy at same point as v (t+dt/2)

  x_temp = nod%x  
  nod%x_prev = nod%x 
  ! if (x_at_midtime  .eqv. .True. ) then
  ! nod%x = nod%x + (1.0d0-gamma)* dt * nod%v
  ! end if 
  !!!!! ACCORDING TO BENSON; STRAIN RATES ARE CALCULATED AT t+1/2dt
  !!!!! ALTERNATED WITH POSITIONS AND FORCES
  call calculate_element_Jacobian()  
  call calc_elem_vol
  call calculate_element_derivMat() !!! WITH NEW SHAPE
  
  call calc_def_grad  
  !call 
  call calc_polar_urmat
  call cal_elem_strain_inc_from_umat !
  
  call calc_elem_pressure_from_strain(domi%mat_K)
  
  !call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  
  !print *, "STRAIN RATE ", elem%str_rate
  !call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES

  nod%x = x_temp

  !!!! SHAPES DERIVATIVES ARE RECALCULATED FOR FORCES CALCULATIONS IN NEW POSITIONS  
  ! call calculate_element_Jacobian()  
  ! call calc_elem_vol
  ! call calculate_element_derivMat() !!! WITH NEW SHAPE  

  
  

  do e=1,elem_count
  if (elem%gausspc(e) > 1) then
      call calculate_element_dhxy0
    end if
  end do
  
  call calc_elem_density
  !print *, "Element density ", elem%rho(:,:)
  !call calc_elem_pressure
  print *, "Element pressure ", elem%pressure(:,:)
  
  
  !call CalcStressStrain(dt)
  call CalcStress (dt)
  
  !print *, "ELEMENT STRESSES ",elem%sigma(:,:,:,:)
  
  !print *, "VELOCITY", nod%v(:,:)  
  call calc_hourglass_forces
  call cal_elem_forces
  call assemble_forces

  !print *, "Element strain rates" 
  do e=1,elem_count
    do gp=1, elem%gausspc(e)
      print *, elem%str_rate(e,gp,:,:)
    end do
  end do

  fext_glob = 0.0d0 !!!ELEMENT 1, node 3,
  
  !print *, "global int forces ", rint_glob(3,:)
  
  
    do n=1,node_count
      do d=1,dim
        nod%a(n,d) =  (fext_glob(n,d)-rint_glob(n,d))/mdiag(n) 
      end do 
    end do
  call impose_bca
  
  nod%a = nod%a - alpha * prev_a
  nod%a = nod%a / (1.0d0 - alpha)
  
  print *, "Accel "
  do n=1,node_count
    print *, nod%a(n,:)
  end do
  
  nod%v = nod%v + gamma * dt * nod%a   
  call impose_bcv !!!REINFORCE VELOCITY BC

  !u = u + beta * nod%v * dt
  u = u + beta * dt * dt * nod%a   
  nod%u = nod%u + u
  nod%x = nod%x + u
  
  !call AverageData(elem%rho(:,1),nod%rho(:))  
  prev_a = nod%a
  time = time + dt
  end do !time ----------------------------------------------------------------------------------------------------------
  
  

  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains
  
end subroutine SolveGreenNag

end module SolverGreenNag