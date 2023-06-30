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

!!!!!!!!!!!!!!!!!!!!!!!!!! Explicit time integration in finite element
! ! method for structural dynamic, wave
! ! propagation and contact-impact problems:
! ! A recent progress
!!!!!! Radek Kolmana , Jos´e Gonz´alezb , J´an Kopaˇckaa , Michal Mraˇckoa
!!!! Evaluate force residual: r
!!!! 0 = fext(t = 0) − Ku0
!!!!!Compute acceleration: u¨0 = M−1 r_0
!!! HERE CORRECT v t - 1/2

module SolverLeapfrog
use ModPrecision, only : fp_kind

contains 

! inline void Particle::Move_Leapfrog(Mat3_t I, double dt)
! {
	! if (FirstStep)
	! {
		! Densitya = Density - dt/2.0*dDensity;
		! va = v - dt/2.0*a;
	! }
	! Densityb = Densitya;
	! Densitya += dt*dDensity;
	! Density = (Densitya+Densityb)/2.0;
	! vb = va;
	! va += dt*a;
	! v = (va + vb)/2.0;
	! x += dt*va;

! inline void Particle::Mat2Leapfrog(double dt)
! {
	! Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	! // Jaumann rate terms
	! Mat3_t RotationRateT,SRT,RS;
	! Trans(RotationRate,RotationRateT);
	! Mult(ShearStress,RotationRateT,SRT);
	! Mult(RotationRate,ShearStress,RS);

	! // Elastic prediction step (ShearStress_e n+1)
	! if (FirstStep)
		! ShearStressa	= -dt/2.0*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;

	! ShearStressb	= ShearStressa;
	! ShearStressa	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressa;


subroutine SolveLeapfrog (tf, dt)
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
 real(fp_kind), dimension(node_count,dim) :: v_med !! v_t+1/2
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
  v_med = 0.0d0
  nod%u(:,:) = 0.0d0
  debug_mode = .true.
  first_step  = .true.
  
  !!!!!!!!!!!!!!! IF EXTERNAL FORCES (AND IF NOT?????, IF BCs ARE ONLY VELOCITY??
  !!!!!!!!!!!!!! CALCULATE Ku0 = RINT0, Initial internal forces
  print *, "Assemblying forces..."
  call assemble_forces()
  do n=1,node_count
      nod%a(n,:) = (fext_glob(n,:)-rint_glob(n,:))/mdiag(n) 
      print *, "fext n ", n, fext_glob(n,:)
  end do
  ! do n=1,node_count
    ! print *, "Initial accel ", n, "a ", nod%a(n,:)  
  ! end do  
  
  nod%v = nod%v - dt * 0.5 * nod%a   !!!!!!!!!!!!!!!!!!v(t -dt/2)
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
  print *,"main loop-------------------------------------"
  do while (time < tf)
    step = step + 1
    print *, "Time: ", time, ", step: ",step, "---------------------------------------------------------"

  
  nod%v = nod%v + dt * nod%a   
  call impose_bcv !!!REINFORCE VELOCITY BC


  !!(3) The velocity is integrated to give the displacement at tn+1.
  nod%u = nod%u +  nod%v * dt
  nod%x = nod%x + nod%u

  !!!! JACOBIAN TO UPDATE SHAPE right after CHANGE POSITIONS
  !!!! IN ORDER TO CALC VOL
  call calculate_element_Jacobian()  
  call calc_elem_vol
  call calculate_element_derivMat() !!! WITH NEW SHAPE

  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES

  ! (7) Based on the density and energy at t_n+l, the pressure is calculated from the equation of
  ! state.
  !!! THIS IS CALCULATE NOW IN ORDER TO UPDATE STRESS WITH CURRENT PRESSURE
  do e=1,elem_count
  if (elem%gausspc(e) > 1) then
      call calculate_element_dhxy0
    end if
  end do
  call calc_elem_density
  call calc_elem_pressure

  call CalcStressStrain(dt)
  
  ! (5) The artificial shock viscosity and hourglass viscosity are calculated from un+1/2. ATTENTION
  call calc_hourglass_forces

    !print *, "det EXT(e,gp)", elem%detJ(:,:)
    do n=1,elem_count
      if (elem%gausspc(n) .eq. 8) then !!!! ELSE IS CONSTANT
        call calculate_element_shapeMat() !AND MASS
      end if
    end do

    call cal_elem_forces()
    call assemble_forces()
    
    !print *, "calc accel "
    do n=1,node_count
      do d=1,dim
        nod%a(n,d) =  (fext_glob(n,d)-rint_glob(n,d))/mdiag(n) 
      end do 
    end do
    call impose_bca
    
    !! IF OUTPUT CORRECT VEL
  time = time + dt
  end do !time
  
  !At the end
  nod%v = nod%v + dt/2. * nod%a   

end subroutine SolveLeapfrog

end module SolverLeapfrog