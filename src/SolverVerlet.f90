module Solver
use ModPrecision, only : fp_kind

contains 

subroutine SolveVerlet (tf, dt)
  use omp_lib
  use Matrices
  implicit none
  integer :: n, d, iglob
  
  real(fp_kind),intent(in)::tf, dt
  
  !NODAL CALCULATION
  
  !Predictor 
  !uest_n+1 = un + dt v_n + dt2/2 a_n
  !Estimate u and vel from previous steps
  !Solve eqns of motion at t_n+1 = tn +dt
  !Calculate a from  M dacc = fext (tn+1) - fint(uest, vest) -fcont
  !Calculate Lumped matrix
  call assemble_mass_matrix()
  !Diagonalize
  !SIMPLEST FORM, ROW SUM 
  iglob = 1
  do while (iglob .le. node_count x dim)
    n = 1
    do while (n .le. node_count x dim) !column
      if (n .neq. iglob)
       !m_glob(iglob,iglob) = m_glob(iglob,iglob)
    end do !col
  iglob = iglob + 1
  end do 
  !Calculate Nodal accelerations a(t+dt) from rext(t)-rint(t)-fcont
  n = 1
  do while (n .le. node_count)
    d = 1
    do while (d .le. 2)
      iglob = (n-1) * dim + d !TODO: CALL THIS IN A FUNCTION
      acc_glob(n,:) = rint_glob(n,d)/m_glob(iglob,iglob) 
      d = d + 1
    end do
    !Solve motion at t_n+1
    !Update vel with CURRENT ACCELERATION
    v_glob(n,:) = v_glob(n,:) + 
    
    n = n + 1
  end do !Node


end subroutine SolveVerlet

end module Solver