module SolverVerlet
use ModPrecision, only : fp_kind

contains 

subroutine SolveVerlet (tf, dt)
  use omp_lib
  use Matrices
  implicit none
  integer :: n, d, iglob
  
  real(fp_kind),intent(in)::tf, dt
  real(fp_kind), dimension(node_count*dim) :: mdiag !!DIAGONALIZATION COULD BE DONE INSIDE ACC CALC  
  real(fp_kind), dimension(dim) :: prev_acc
  
  
  call calculate_element_matrices()!ATTENTION: THIS CALCULATES KNL AND KL AND THIS REQUIRES UPDATE CAUCHY STRESS TAU
  !NODAL CALCULATION
  
  !Predictor 
  !uest_n+1 = un + dt v_n + dt2/2 a_n
  !Estimate u and vel from previous steps
  !Solve eqns of motion at t_n+1 = tn +dt
  !Calculate a from  M dacc = fext (tn+1) - fint(uest, vest) -fcont
  !Calculate Lumped matrix
  call assemble_mass_matrix()
  call assemble_int_forces()
  ! !Diagonalize
  ! !SIMPLEST FORM, ROW SUM 
  !print *, "Calc mdiag"
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
  !Calculate positions with PREVIOUS ACCELERATIONS
  n = 1
  print *, "calculating positions "
  do while (n .le. node_count)
    !print *, "n ", n, "acc ", nod%v(n,:)
    !nod%x(n,:) = nod%x(n,:) + nod%v(n,:) * dt + nod%a(n,:) * dt * dt    
    !nod%x(n,:) = nod%x(n,:) +  nod%v(n,:) * dt + nod%a(n,:) * dt * dt    
    n = n + 1
  end do !Node    
  !Calculate Nodal accelerations a(t+dt) from rext(t)-rint(t)-fcont
  print *, "calculating accel"
  n = 1
  do while (n .le. node_count)
    !print *, "node ", n 
    prev_acc(:) = nod%a(n,:)
    d = 1
    do while (d .le. 2)
      !print *, "dim ", d
      iglob = (n-1) * dim + d !TODO: CALL THIS IN A FUNCTION
      !print *, "iglob ", iglob, "mdiag ", mdiag(iglob)
      
      !nod%a(n,d) = rint_glob(n,d)/mdiag(iglob) 
      d = d + 1
    end do
    !Solve motion at t_n+1
    !Update vel with CURRENT ACCELERATION
    !nod%v(n,:) = nod%v(n,:) + dt * 0.5 * (nod%a(n,:) + prev_acc(:)) !THIS CRASH
    
    n = n + 1
  end do !Node


end subroutine SolveVerlet

end module SolverVerlet