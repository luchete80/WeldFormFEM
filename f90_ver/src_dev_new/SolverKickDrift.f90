! ////////////////////////////////////////////////////////
! ////////////////////// BASED ON RANDLES AND LIBERSKY(1996):
! /////////////////////////////////////////////////////////
! /// Randles and Libersky calculate density from current velocity, here is from the velocity at t+1/2
    ! // // // 1 CalcAccel(); //Nor density or neither strain rates
    ! // // // if (nonlock_sum)AccelReduction();

    ! // // // if (contact) CalcContactForcesWang();

    ! // // // 3. Particles[i]->v += Particles[i]->a*deltat/2.*factor;

    ! // // // 4. //If density is calculated AFTER displacements, it fails
    ! // // // CalcDensInc(); //TODO: USE SAME KERNEL?
    ! // // // Particles[i]->Density += deltat*Particles[i]->dDensity*factor;
    ! // // // 5. x += (Particles[i]->v + Particles[i]->VXSPH)*deltat*factor;
    
    ! // // // 6. Particles[i]->v += Particles[i]->a*deltat/2.*factor;
    ! // // // 7. CalcRateTensors();  //With v and xn+1
    ! // // // 8. Particles[i]->CalcStressStrain(deltat); //Uses density  


module SolverKickDrift
use ModPrecision, only : fp_kind

contains 

subroutine SolveKickDrift (tf, dt)
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
  
  
  do n=1,node_count
    print *, "Initial accel ", n, "a ", nod%a(n,:)  
  end do  
  
  !!!! IF ONLY ARE SET bcv 
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
  print *,"------------------------------------------------------------------------------------------------"
  print *,"main loop, VERLET -------------------------------------------------------------------------------"
  do while (time < tf)
    step = step + 1
    print *, "Time: ", time, ", step: ",step, "---------------------------------------------------------"

  !!!!!!!!!!!!!!! IF EXTERNAL FORCES (AND IF NOT?????, IF BCs ARE ONLY VELOCITY??
  !!!!!!!!!!!!!! CALCULATE Ku0 = RINT0, Initial internal forces
  print *, "Assemblying forces..."
  fext_glob = 0.0d0 !!!ELEMENT 1, node 3,
  call cal_elem_forces
  call assemble_forces()
  do n=1,node_count
      nod%a(n,:) = (fext_glob(n,:)-rint_glob(n,:))/mdiag(n) 
      print *, "fext n ", n, fext_glob(n,:)
  end do
  call impose_bca
  
  
  ! if (time < 100.0d0*dt) then
    ! nod%bcv(5:8,3) = -0.1 * time/(10.0d0*dt)
    ! !nod%bcv(3:4,2) = -0.1 * time/(100.0d0*dt)
  ! else 
    ! nod%bcv(5:8,3) = -0.1
    ! !nod%bcv(3:4,2) = -0.1
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

  call calc_elem_density

  ! !!(3) The velocity is integrated to give the displacement at tn+1.
  nod%u = nod%u +  nod%v * dt!/2.0  
  nod%x = nod%x + nod%u

    ! !!!!! ACCORDING TO BENSON; STRAIN RATES ARE CALCULATED AT t+1/2dt
  ! !!!!! ALTERNATED WITH POSITIONS AND FORCES
  call calculate_element_Jacobian()  
  call calc_elem_vol
  call calculate_element_derivMat() !!! WITH NEW SHAPE
  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES
  
  nod%v = nod%v + dt/2.0 * nod%a   
  call impose_bcv !!!REINFORCE VELOCITY BC


    ! ! !!!!! ACCORDING TO BENSON; STRAIN RATES ARE CALCULATED AT t+1/2dt
  ! ! !!!!! ALTERNATED WITH POSITIONS AND FORCES
  ! call calculate_element_Jacobian()  
  ! call calc_elem_vol
  ! call calculate_element_derivMat() !!! WITH NEW SHAPE
  ! call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  ! call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES



  !!!! SHAPES DERIVATIVES ARE RECALCULATED FOR FORCES CALCULATIONS IN NEW POSITIONS  
  ! call calculate_element_Jacobian()  
  ! call calc_elem_vol
  ! call calculate_element_derivMat() !!! WITH NEW SHAPE  

  
  

  ! do e=1,elem_count
  ! if (elem%gausspc(e) > 1) then
      ! call calculate_element_dhxy0
    ! end if
  ! end do
  


  call calc_elem_pressure
  call CalcStressStrain(dt)


  time = time + dt
  end do !time

  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains
  
end subroutine SolveKickDrift

end module SolverKickDrift