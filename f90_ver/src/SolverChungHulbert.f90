
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

module SolverChungHulbert
use ModPrecision, only : fp_kind

contains 

!!!COMPUTE STEP SIZE
!!!maximumFrequency = 2.0 / model->computeCourantTimeStep();
!!!timeStep = _timeStepSafetyFactor * _omegaS / maximumFrequency;


subroutine SolveChungHulbert (domi, tf, dt)
  use omp_lib
  use Matrices
  use Mechanical
  
  implicit none
  integer :: n, d, iglob, step, e, gp
	integer :: last_out
  integer :: step_out
  
  logical :: first_step, x_at_midtime
  logical :: debug_mode 
  real(fp_kind),intent(in)::tf, dt
  type (dom_type), intent (in) :: domi
  
  real(fp_kind), dimension(node_count) :: mdiag !!DIAGONALIZATION COULD BE DONE INSIDE ACC CALC  
  real(fp_kind), dimension(dim) :: prev_acc
  real(fp_kind), dimension(node_count,dim) :: u, prev_a, x_temp
  
  real(fp_kind) :: alpha, beta, gamma, rho_b , omega!!! CHUNG HULBERT PARAMETERS
 
  real(fp_kind), dimension(nodxelem,dim) :: xtest

  call set_edof_from_elnod()
  
  !! IN AXISYMM CASE: calc matm which depends on Jacobian, jacobian BUT JACOBIAN NEEDS RADIUS
  call calculate_element_Jacobian  
  print *, "shape mat"
  if (bind_dom_type .eq. 3) then
    call calculate_element_shapeMat() !ONLY FOR VOLUMETRIC CALCS
    print *, "shape mat", elem%math(1,1,:,:)
    print *, "calc radis"
    call calculate_element_radius()
    call calculate_element_MassMat ()
  end if
  
  call calc_elem_vol !!!! In order to define initial volume
  call calculate_element_derivMat()
  elem%vol_0(:) = elem%vol(:)
  !print *,"Element Initial Vol"
  ! do n = 1, elem_count
    ! print *, elem%vol(n)
  ! end do    
  call calc_nodal_masses()
    
  ! call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU

  if (bind_dom_type .eq. 3) then !!! AXISYMM
    print *, "Assemblying mass matrix" 
    call assemble_mass_matrix() !!! mglob
  end if 

  ! mdiag(:)=0.0d0

  print *, "Nodal masses"
  print *,  nod%m(:)
  
  calc_m = .False.
  ! print *, "M Diag with mass mat", mdiag
  PRINT *, "Tot Mass: ", tot_mass
  !print *, "Tot mass from mdiag", 
  !!!! ONLY FOR TESTING
  if (bind_dom_type .ne. 3) then
    do n=1, node_count  !column
      mdiag(n) = tot_mass/node_count 
      nod%m(iglob) = mdiag(iglob)
    end do
  else
    if (axisymm_vol_weight .eqv. .true.) then 
      mdiag(:)=0.0d0
      do iglob =1, node_count
        do n=1, node_count  !column
           mdiag(iglob) = mdiag(iglob) + m_glob(iglob,n)
        end do !col
          if (axisymm_vol_weight .neqv. .true.) then 
            nod%m(iglob) = mdiag(iglob)
          end if 
        ! print *,  mdiag(iglob)
      end do 
    else  !! AXISYMM AREA WEIGHT
    do n=1, node_count  !column
      !!According to goudreau (19), Benson 
      mdiag(n) = nod%m(n)
      ! nod%m(iglob) = mdiag(iglob)
    end do    
    end if
  end if

  ! print *, "M Diag with node avg", mdiag
    print *, "Nodal masses"
  print *,  mdiag(:)
  
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
  beta  = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b))
  gamma = 1.5 - alpha;

  print *, "alpha ", alpha
  print *, "beta ", beta
  print *, "gamma ", gamma
  
  !In central difference (leapfrog)
  
  
  elem%shear_stress = 0.0d0 
  time = 0.0  
  step = 0
  
  x_at_midtime = .False.
  
  step_out = 1 !FREQUENCY
	last_out = 0

  print *, "mdiag ", mdiag(n) 
  !! ONLY TO TEST THINGS
  
  print *,"------------------------------------------------------------------------------------------------"
  print *,"main loop, CHUNG HULBERT -----------------------------------------------------------------------"
  do while (time < tf)
    step = step + 1
		if (step > last_out) then 
			print *, "Time: ", time, ", step: ",step, "---------------------------------------------------------"
			last_out = last_out + step_out
		end if 
  ! if (time < 100.0d0*dt) then
    ! nod%bcv(5:8,3) = -1.0 * time/(10.0d0*dt)
    ! !nod%bcv(3:4,2) = -0.1 * time/(100.0d0*dt)
  ! else 
    ! nod%bcv(5:8,3) = -1.0
    ! !nod%bcv(3:4,2) = -0.1
  ! end if 
  
    ! do n=1,elem_count
      ! if (elem%gausspc(n) .eq. 8) then !!!! ELSE IS CONSTANT
        ! call calculate_element_shapeMat() !AND MASS
      ! end if
    ! end do
  !!! PREDICTION PHASE
  u = dt * (nod%v + (0.5d0 - beta) * dt * prev_a)
  !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
  nod%v = nod%v + (1.0d0-gamma)* dt * prev_a
  nod%a = 0.0d0
  
  call impose_bcv !!!REINFORCE VELOCITY BC
  !print *, "veloc", nod%v 
  ! nod%u = nod%u +  nod%v * dt!/2.0  
  ! nod%x = nod%x + nod%u             !! EVALUATE dHdxy at same point as v (t+dt/2)

  x_temp = nod%x  
  if (x_at_midtime  .eqv. .True. ) then
  nod%x = nod%x + (1.0d0-gamma)* dt * nod%v
  end if 
  !!!!! ACCORDING TO BENSON; STRAIN RATES ARE CALCULATED AT t+1/2dt
  !!!!! ALTERNATED WITH POSITIONS AND FORCES
  call calculate_element_Jacobian()  
  call calc_elem_vol
  call calculate_element_derivMat() !!! WITH NEW SHAPE
  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES

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
  !!!call calc_elem_pressure
  
  call cal_elem_strain_inc_from_str_rate (dt)
  call calc_elem_pressure_from_strain(domi%mat_K)  
  
  !print *, "Element pressure ", elem%pressure(:,:)

	!! TODO: MODIFY THIS
	! if (dim .eq. 3) then 
		 call CalcStressStrain(dt)
  ! else
		! call Calc_Elastic_Stress(domi, dt) !!!ELASTIC_TEST
  ! end if
	!print *, "VELOCITY", nod%v(:,:)  
  call calc_hourglass_forces
  call cal_elem_forces
  call assemble_forces

  !print *, "Element strain rates" 
  ! do e=1,elem_count
    ! do gp=1, elem%gausspc(e)
      ! print *, elem%str_rate(e,gp,:,:)
    ! end do
  ! end do

  fext_glob = 0.0d0 !!!ELEMENT 1, node 3,
  
  !print *, "global int forces ", rint_glob(3,:)
  
  !!!! IF AREA WEIGHTED UPDATE ELEMENT MASS WITH RADIUS
  ! if (bind_dom_type .eq. 3 .and. axisymm_vol_weight .eqv. .false.) then
    ! do n=1,node_count
      
    ! end do
  ! end if

	!$omp parallel do num_threads(Nproc) private (n)  
	do n=1,node_count
		! do d=1,dim
			! nod%a(n,d) =  (fext_glob(n,d)-rint_glob(n,d))/mdiag(n) 
		! end do 
		nod%a(n,:) =  (fext_glob(n,:)-rint_glob(n,:))/mdiag(n) 
	end do
	!$omp end parallel do
	
    
  call impose_bca
  
	!$omp parallel do num_threads(Nproc) private (n)
  do n=1,node_count
		nod%a(n,:) = nod%a(n,:) - alpha * prev_a(n,:)
		nod%a(n,:) = nod%a(n,:) / (1.0d0 - alpha)
		nod%v(n,:) = nod%v(n,:) + gamma * dt * nod%a (n,:)  
	end do
	!$omp end parallel do
  

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

  call CalcEquivalentStress()
  call AverageData(elem%rho(:,1),nod%rho(:))
  call AverageData(elem%sigma_eq(:,1),nod%sigma_eq(:))
  call AverageData(elem%pl_strain(:,1),nod%pl_strain(:))
  
  ! call WriteMeshVTU('output.vtu')
  
end subroutine SolveChungHulbert

end module SolverChungHulbert
