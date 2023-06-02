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
  integer :: n, d, iglob
  
  logical :: first_step 
  
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
  
  first_step  = .true.
  
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
    
  call impose_bca
  if (first_step) then 
    nod%v(n,:) = nod%v(n,:) - dt * 0.5 * nod%a(n,:)
  end if

  !!!!! THIS IS NOT SOLVED AS A COMPLETED STEP (REDUCED VERLET=
  ! (2) The acceleration is integrated to give the velocity at tn+l/2.
  ! !Update vel with CURRENT ACCELERATION
  ! THIS WOULD BE AT ONE STEP
  ! nod%v(n,:) = nod%v(n,:) + dt * 0.5 * (nod%a(n,:) + prev_acc(:)) 
  ! print *,"node vel ", nod%v(n,:)  
  nod%v(n,:) = nod%v(n,:) + dt * 0.5 * nod%a(n,:)
  call impose_bcv !!!REINFORCE VELOCITY BC

  !!(3) The velocity is integrated to give the displacement at tn+1.
  nod%x(n,:) = nod%x(n,:) +  nod%v(n,:) * dt
  
  call calc_elem_vol
  call disassemble_uvele     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
  call cal_elem_strains      !!!!!STRAIN AND STRAIN RATES
  
  ! (4) The constitutive model for the strength of the material is integrated from t to t_n+1 now
  ! that the motion of the material is known.
  
  ! (5) The artificial shock viscosity and hourglass viscosity are calculated from un+1/2. ATTENTION
  call calc_hourglass_forces
  

! (6) The internal energy is updated based on the work done between tn and t_n+1.

! (7) Based on the density and energy at t_n+l, the pressure is calculated from the equation of
! state.
  call calc_elem_density
  call calc_elem_pressure
  
  
  nod%v(n,:) = nod%v(n,:) + dt * 0.5 * nod%a(n,:)
  call impose_bcv !!!REINFORCE VELOCITY BC

  time = time + dt
  end do !time

end subroutine SolveLeapfrog

end module SolverLeapfrog