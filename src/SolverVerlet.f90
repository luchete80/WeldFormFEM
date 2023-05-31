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
  calc_m = .True.
 
  call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU
  call assemble_mass_matrix()
    mdiag(:)=0.0d0
    iglob = 1
    do while (iglob .le. node_count * dim)
      n = 1
      do while (n .le. node_count * dim) !column
         mdiag(iglob) = mdiag(iglob) + m_glob(iglob,n)
         n = n+ 1
      end do !col
    iglob = iglob + 1
    end do   
  calc_m = .False.
  print *, "M Diag ", mdiag
  print *, "m glob", m_glob
  
  nod%u(:,:) = 0.0d0
  
  time = 0.0  
  print *,"main loop"
  do while (time .le. tf)

    call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU
    !NODAL CALCULATION
    
    !Predictor 
    !uest_n+1 = un + dt v_n + dt2/2 a_n
    !Estimate u and vel from previous steps
    !Solve eqns of motion at t_n+1 = tn +dt
    !Calculate a from  M dacc = fext (tn+1) - fint(uest, vest) -fcont
    !Calculate Lumped matrix
    
    call calc_elem_vol()
    call disassemble_uele()     !BEFORE CALLING UINTERNAL AND STRAIN STRESS CALC
    call assemble_int_forces()
    ! !Diagonalize
    ! !SIMPLEST FORM, ROW SUM 
    !print *, "Calc mdiag"
    ! if (calc_m .eqv. .True.) then 
      ! call assemble_mass_matrix()
      ! mdiag(:)=0.0d0
      ! iglob = 1
      ! do while (iglob .le. node_count * dim)
        ! n = 1
        ! do while (n .le. node_count * dim) !column
           ! mdiag(iglob) = mdiag(iglob) + m_glob(iglob,n)
           ! n = n+ 1
        ! end do !col
      ! iglob = iglob + 1
      ! end do 
    ! end if
    
    print *, " mglob ," ,m_glob
    print *, " mdiag ," ,mdiag
    
    call impose_bcv
    !Calculate positions with PREVIOUS ACCELERATIONS
    n = 1
    print *, "calculating positions "
    do while (n .le. node_count)
      !print *, "n ", n, "acc ", nod%v(n,:)
      nod%u(n,:) = nod%u(n,:) + nod%v(n,:) * dt + nod%a(n,:) * dt * dt  
      nod%x(n,:) = nod%x(n,:) + nod%u(n,:)    
      !nod%x(n,:) = nod%x(n,:) +  nod%v(n,:) * dt + nod%a(n,:) * dt * dt    
      print *,"node u ", nod%u(n,:)
      n = n + 1
    end do !Node 
    
    !Again to calculate strains
    call disassemble_uele()
    print *, "u ele 2", elem%uele(2,:,:)
    call cal_elem_strains ()
    
    print *,"strains", elem%strain(:,:,:,:)
    !Calculate Nodal accelerations a(t+dt) from rext(t)-rint(t)-fcont
    print *, "calculating accel"
    
    do n=1,node_count
      !print *, "node ", n 
      prev_acc(:) = nod%a(n,:)
      d = 1
      do while (d .le. 2)
        !print *, "dim ", d
        iglob = (n-1) * dim + d !TODO: CALL THIS IN A FUNCTION
        print *, "iglob ", iglob, "mdiag ", mdiag(iglob)
        
        nod%a(n,d) = rint_glob(n,d)/mdiag(iglob) 
        print *,"node acc ", nod%a(n,:), "rint ", rint_glob(n,:)
        d = d + 1
      end do
      !Solve motion at t_n+1
      !Update vel with CURRENT ACCELERATION
      nod%v(n,:) = nod%v(n,:) + dt * 0.5 * (nod%a(n,:) + prev_acc(:)) 
      print *,"node vel ", nod%v(n,:)

    end do !Node
    call impose_bca
    
    !REINFORCE bc velocity
    !TODO: SEPARATE FUNCTION
    call impose_bcv
    
    !calc_m = .False.
  time = time + dt
end do !time

end subroutine SolveVerlet

end module SolverVerlet