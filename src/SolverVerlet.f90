!!!!BENSON 1992
    !!calculate v)t+1/2)

    ! CalcAccel(); //Nor density or neither strain rates
    ! AccelReduction();
    ! GeneralAfter(*this); //Fix free accel

  ! if (contact) CalcContactForcesWang();
    
    ! Particles[i]->v += Particles[i]->a*deltat/2.*factor;
    ! MoveGhost();   
    ! GeneralAfter(*this);//Reinforce BC vel   

    ! CalcDensInc(); //TODO: USE SAME KERNEL?
      ! Particles[i]->Density += deltat*Particles[i]->dDensity*factor;
      
      ! du = (Particles[i]->v + Particles[i]->VXSPH)*deltat*factor;
      ! Particles[i]->Displacement += du;
      ! Particles[i]->x += du;

    ! for (size_t i=0; i<Particles.Size(); i++){
      ! Particles[i]->v += Particles[i]->a*deltat/2.*factor;

    ! GeneralAfter(*this);
    ! CalcRateTensors();  //With v and xn+1
    ! Particles[i]->CalcStressStrain(deltat); //Uses density  

		! clock_beg = clock();        
    ! CalcKinEnergyEqn();    
    ! CalcIntEnergyEqn();    
    ! UpdateContactParticles(); //Updates normal and velocities
		

module SolverVerlet
use ModPrecision, only : fp_kind

contains 

subroutine SolveVerlet (tf, dt)
  use omp_lib
  use Matrices
  use Mechanical
  
  implicit none
  integer :: n, d, iglob, step, e, gp
  
  logical :: first_step 
  logical :: debug_mode 
  real(fp_kind),intent(in)::tf, dt
  
  real(fp_kind), dimension(node_count) :: mdiag !!DIAGONALIZATION COULD BE DONE INSIDE ACC CALC  
  real(fp_kind), dimension(dim) :: prev_acc
 
  real(fp_kind), dimension(nodxelem,dim) :: xtest

  call set_edof_from_elnod()
  
  call calculate_element_Jacobian()
  call calculate_element_shapeMat() !AND MASS
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
  print *, "Assemblying forces..."
  call assemble_forces()
  do n=1,node_count
      nod%a(n,:) = (fext_glob(n,:)-rint_glob(n,:))/mdiag(n) 
      print *, "fext n ", n, fext_glob(n,:)
  end do
  call impose_bca
  
  ! do n=1,node_count
    ! print *, "Initial accel ", n, "a ", nod%a(n,:)  
  ! end do  
  
  !nod%v = nod%v - dt * 0.5 * nod%a   !!!!!!!!!!!!!!!!!!v(t -dt/2)
  call impose_bcv

  ! do n=1,node_count
    ! print *, "Initial v nod ", n, "v ", nod%v(n,:)  
  ! end do  
  
  !!!! IS THERE ANY STRESS?
  elem%sigma (:,:,:,:) = 0.0d0 !!!! FOR INT FORCES (elem%f_int(e,gp,d,d)) CALCULATION
  !elem%f_int (:,:,:)   = 0.0d0 !!!! I Ncal_elem_forces

  elem%shear_stress = 0.0d0 
  time = 0.0  
  step = 0
  print *,"main loop---------------------------------------------------------------------------------------"
  do while (time < tf)
    step = step + 1
    print *, "Time: ", time, ", step: ",step, "---------------------------------------------------------"

    !print *, "det EXT(e,gp)", elem%detJ(:,:)
    do n=1,elem_count
      if (elem%gausspc(n) .eq. 8) then !!!! ELSE IS CONSTANT
        call calculate_element_shapeMat() !AND MASS
      end if
    end do
    
    
    !!! TEST 
    !xtest(:,:) = []
    
    ! call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU
    ! !NODAL CALCULATION

    ! call assemble_int_forces()
  ! (1) Knowing the stress, pressure, hourglass forces and shock viscosity at t” in each zone or
  ! element, the forces at the nodes are calculated. The accelerations of the nodes are
  ! calculated by dividing the nodal forces by the nodal masses.   
    !print *, "calc elem forces "
    call cal_elem_forces()

    if (debug_mode .eqv. .true.) then    
    do e=1, elem_count
      do gp = 1, elem%gausspc(e)
      !print *, "elem%dHxy_detJ(e,gp,1", elem%dHxy_detJ(e,gp,1,:)
      !print *, "elem%dHxy_detJ(e,gp,2", elem%dHxy_detJ(e,gp,2,:)
        do n=1, nodxelem
          print *, "Node ", n, "Force ", elem%f_int(e,n,:) 
        end do! nod x elem
      end do !gp
    end do!elem      
    end if  
    !print *, "assemble int forces "
    call assemble_forces()
    
    !print *, "calc accel "
    do n=1,node_count
      do d=1,dim
        nod%a(n,d) =  (fext_glob(n,d)-rint_glob(n,d))/mdiag(n) 
      end do 
    end do
  if (debug_mode .eqv. .true.) then
    do n=1,node_count
      print *, "glob res forces ", (fext_glob(n,:)-rint_glob(n,:))
    end do 
  end if
  call impose_bca
  
  !!!!! THIS IS NOT SOLVED AS A COMPLETED STEP (REDUCED VERLET=
  ! (2) The acceleration is integrated to give the velocity at tn+l/2.
  ! !Update vel with CURRENT ACCELERATION
  ! THIS WOULD BE AT ONE STEP
  ! nod%v(n,:) = nod%v(n,:) + dt * 0.5 * (nod%a(n,:) + prev_acc(:)) 
  ! print *,"node vel ", nod%v(n,:)  
  nod%v = nod%v + dt * nod%a   
  call impose_bcv !!!REINFORCE VELOCITY BC

  if (debug_mode .eqv. .true.) then
  do n=1,node_count
    print *, "nod ", n, "a ", nod%a(n,:)  
  end do  
  end if 
  !!(3) The velocity is integrated to give the displacement at tn+1.
  nod%u = nod%u +  nod%v * dt
  nod%x = nod%x + nod%u

  if (debug_mode .eqv. .true.) then
  do n=1,node_count
    print *, "nod ", n, "v ", nod%v(n,:)  
  end do  
  print *, "delta t", dt
  do n=1,node_count
    print *, "nod ", n, "x ", nod%x(n,:)  
  end do  
  end if
  !!!! JACOBIAN TO UPDATE SHAPE right after CHANGE POSITIONS
  !!!! IN ORDER TO CALC VOL
  call calculate_element_Jacobian()  
  call calc_elem_vol
  call calculate_element_derivMat() !!! WITH NEW SHAPE
  
  if (debug_mode .eqv. .true.) then
    print *,"Element Vol"
    do n = 1, elem_count
      print *, elem%vol(n)
    end do
    print *,"Element Mass"
    do n = 1, elem_count
      print *, elem%mass(n)
    end do
  end if
  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES

  ! (7) Based on the density and energy at t_n+l, the pressure is calculated from the equation of
  ! state.
  !!! THIS IS CALCULATE NOW IN ORDER TO UPDATE STRESS WITH CURRENT PRESSURE
  call calc_elem_density
  call calc_elem_pressure
  
  if (debug_mode .eqv. .true.) then
    print *,"Element Density"
    do n = 1, elem_count
      print *, elem%rho(n,:)
    end do
    
    print *,"Element pressure"
    do n = 1, elem_count
      print *, elem%pressure(n,:)
    end do
  
  end if 
  
  ! (4) The constitutive model for the strength of the material is integrated from t to t_n+1 now
  ! that the motion of the material is known.
  print *, "Calc stresses "
  call CalcStressStrain(dt)
  
  ! (5) The artificial shock viscosity and hourglass viscosity are calculated from un+1/2. ATTENTION
  call calc_hourglass_forces
  if (debug_mode .eqv. .true.) then
    do n=1,elem_count
    print *, "hourglass forces ", elem%hourg_nodf(n,:,:)
    end do
  end if
! (6) The internal energy is updated based on the work done between tn and t_n+1.

  call AverageData(elem%rho(:,1),nod%rho(:))
  !if (debug_mode .eqv. .true.) then
    do n=1,node_count
      print *, "nod ", n, "Disp ", nod%u(n,:)  
    end do  
  !end if
  ! print *, "nod v", nod%v(:,:)
  

  time = time + dt
  end do !time

end subroutine SolveLeapfrog

end module SolverLeapfrog