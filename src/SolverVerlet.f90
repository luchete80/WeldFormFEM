module Solver
use ModPrecision, only : fp_kind

contains 

subroutine SolveVerlet (tf, dt)
  use omp_lib
  use Matrices
  implicit none
  integer :: n, d, iglob
  
  real(fp_kind),intent(in)::tf, dt
  
  real(fp_kind), allocatable, dimension(:) :: prev_acc
  
  allocate (prev_acc(dim))
  
  !NODAL CALCULATION
  
  !Predictor 
  !uest_n+1 = un + dt v_n + dt2/2 a_n
  !Estimate u and vel from previous steps
  !Solve eqns of motion at t_n+1 = tn +dt
  !Calculate a from  M dacc = fext (tn+1) - fint(uest, vest) -fcont
  !Calculate Lumped matrix
  call assemble_mass_matrix()
  call assemble_int_forces()
  !Diagonalize
  !SIMPLEST FORM, ROW SUM 
  iglob = 1
  do while (iglob .le. node_count * dim)
    n = 1
    do while (n .le. node_count * dim) !column
      if (n .ne. iglob) then
       !m_glob(iglob,iglob) = m_glob(iglob,iglob)
      end if
    end do !col
  iglob = iglob + 1
  end do 
  !Calculate positions with PREVIOUS ACCELERATIONS
  n = 1
  do while (n .le. node_count)
    nod%x(n,:) = nod%x(n,:) + nod%v(n,:) * dt + 0.5d0 * nod%a(n,:) * dt    
    n = n + 1
  end do !Node    
  !Calculate Nodal accelerations a(t+dt) from rext(t)-rint(t)-fcont
  n = 1
  do while (n .le. node_count)
    d = 1
    do while (d .le. 2)
      iglob = (n-1) * dim + d !TODO: CALL THIS IN A FUNCTION
      prev_acc = nod%a(n,:)
      nod%a(n,:) = rint_glob(n,d)/m_glob(iglob,iglob) 
      d = d + 1
    end do
    !Solve motion at t_n+1
    !Update vel with CURRENT ACCELERATION
    nod%v(n,:) = nod%v(n,:) + 0.5d0 * (nod%a(n,:) + prev_acc)*dt
    
    n = n + 1
  end do !Node


end subroutine SolveVerlet

end module Solver