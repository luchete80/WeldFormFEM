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

subroutine SolveVerlet (domi, tf, dt) !!!! TODO: REPLACE DOMI FOR MATERIAL
  use omp_lib
  use Matrices
  use Mechanical
  
  implicit none
  integer :: n, d, iglob, step, e, gp
  
  logical :: first_step 
  logical :: debug_mode 
  real(fp_kind),intent(in)::tf, dt
  type (dom_type), intent (in) :: domi
  
  real(fp_kind), dimension(node_count) :: mdiag !!DIAGONALIZATION COULD BE DONE INSIDE ACC CALC  
  real(fp_kind), dimension(dim) :: prev_acc
 
  real(fp_kind), dimension(nodxelem,dim) :: xtest

  call set_edof_from_elnod()
  
  call calculate_element_Jacobian
  print *, "shape mat"
  !call calculate_element_shapeMat() !AND MASS
  call calc_elem_vol !!!! In order to define initial volume
  call calculate_element_derivMat()
  elem%vol_0(:) = elem%vol(:)
  print *,"Element Initial Vol"
  do n = 1, elem_count
    print *, elem%vol(n)
  end do    
    
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

  !!!! IS THERE ANY STRESS?
  elem%sigma (:,:,:,:) = 0.0d0 !!!! FOR INT FORCES (elem%f_int(e,gp,d,d)) CALCULATION
  elem%f_int (:,:,:)   = 0.0d0 !!!! I Ncal_elem_forces
  
  
  elem%sigma = 0.0d0
  fext_glob = 0.0d0
  !!!!!!!!!!!!!!! IF EXTERNAL FORCES (AND IF NOT?????, IF BCs ARE ONLY VELOCITY??
  !!!!!!!!!!!!!! CALCULATE Ku0 = RINT0, Initial internal forces
  print *, "Assemblying forces..."
  call cal_elem_forces
  call assemble_forces()
  do n=1,node_count
      nod%a(n,:) = (fext_glob(n,:)-rint_glob(n,:))/mdiag(n) 
      print *, "fext n ", n, fext_glob(n,:)
      print *, "rint_glob", rint_glob(n,:)
  end do
  call impose_bca
  
  do n=1,node_count
    print *, "Initial accel ", n, "a ", nod%a(n,:)  
  end do  
  
  !!!! IF ONLY ARE SET bcv 
  call impose_bcv

  ! do n=1,node_count
    ! print *, "Initial v nod ", n, "v ", nod%v(n,:)  
  ! end do  
  


  elem%shear_stress = 0.0d0 
  time = 0.0  
  step = 0
  print *,"------------------------------------------------------------------------------------------------"
  print *,"main loop, VERLET -------------------------------------------------------------------------------"
  do while (time < tf)
    step = step + 1
    print *, "Time: ", time, ", step: ",step, "---------------------------------------------------------"

  ! if (time < 5.0e-3) then
    ! nod%bcv(5:8,3) = 0.1 * time/(10.0d0*dt)
    ! !nod%bcv(3:4,2) = 0.1 * time/(100.0d0*dt)
  ! else 
    ! nod%bcv(5:8,3) = 0.1
    ! !nod%bcv(3:4,2) = 0.1
  ! end if 
  
    do n=1,elem_count
      if (elem%gausspc(n) .eq. 8) then !!!! ELSE IS CONSTANT
        call calculate_element_shapeMat() !AND MASS
      end if
    end do
    
  !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
  nod%v = nod%v + dt/2.0 * nod%a   
  call impose_bcv !!!REINFORCE VELOCITY BC
  print *, "veloc", nod%v 
  ! nod%u = nod%u +  nod%v * dt!/2.0  
  ! nod%x = nod%x + nod%u             !! EVALUATE dHdxy at same point as v (t+dt/2)

  
  !!!!! ACCORDING TO BENSON; STRAIN RATES ARE CALCULATED AT t+1/2dt
  !!!!! ALTERNATED WITH POSITIONS AND FORCES
  nod%u = nod%u +  nod%v * dt/2.0  
  nod%x = nod%x + nod%u
  call calculate_element_Jacobian()  
  call calc_elem_vol
  call calculate_element_derivMat() !!! WITH NEW SHAPE
  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES

  ! !!(3) The velocity is integrated to give the displacement at tn+1.
  nod%u = nod%u +  nod%v * dt/2.0  
  nod%x = nod%x + nod%u

  !!!! SHAPES DERIVATIVES ARE RECALCULATED FOR FORCES CALCULATIONS IN NEW POSITIONS  
  call calculate_element_Jacobian()  
  call calc_elem_vol
  call calculate_element_derivMat() !!! WITH NEW SHAPE  

  
  

  do e=1,elem_count
  if (elem%gausspc(e) > 1) then
      call calculate_element_dhxy0
    end if
  end do
  
  call calc_elem_density
  print *, "Element density ", elem%rho(:,:)
  call calc_elem_pressure
  print *, "Element pressure ", elem%pressure(:,:)

  call CalcStressStrain(dt)
  ! print *, "VELOCITY", nod%v(:,:)  
  elem%hourg_nodf(:,:,:) = 0.0d0
  call calc_hourglass_forces
  
  ! !!!! THIS IS FOR SHOCK
  ! print *, "domi%mat_K", domi%mat_K
  call calc_elem_wave_speed(domi%mat_K)
  call calc_elem_shock_visc(dt)
  !!!!
  
  call cal_elem_forces
  call assemble_forces

  print *, "Element strain rates" 
  do e=1,elem_count
    do gp=1, elem%gausspc(e)
      print *, elem%str_rate(e,gp,:,:)
    end do
  end do

  fext_glob = 0.0d0 !!!ELEMENT 1, node 3,

  
  print *, "global int forces ", rint_glob(3,:)
  
    do n=1,node_count
      do d=1,dim
        nod%a(n,d) =  (fext_glob(n,d)-rint_glob(n,d))/mdiag(n) 
      end do 
    end do
  call impose_bca
  
  nod%v = nod%v + dt/2.0 * nod%a   
  call impose_bcv !!!REINFORCE VELOCITY BC

  
  !call AverageData(elem%rho(:,1),nod%rho(:))  

  time = time + dt
  end do !time ----------------------------------------------------------------------------------
  

  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains
  
end subroutine SolveVerlet

end module SolverVerlet