module Mechanical
use Domain
use ElementData
use NodeData

contains
subroutine cal_elem_strains ()
  integer :: e
  integer :: i,j,k
  real(fp_kind), dimension(2) :: r, s
  
 
  r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
  s(1) = r(1)          ; s(2) =  r(2)
  e = 1
  do while (e <= elem_count)  
    !Is only linear matrix?    
    elem%strain(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%uele (e,:,:)) 

    e = e + 1
  end do
end subroutine

subroutine impose_bcv
  implicit none
  integer :: n, d
  n = 1
  do while (n .le. node_count)    
    d = 1
    do while (d .le. 2)
      if (nod%is_bcv(n,d) .eqv. .true.) then
        nod%v(n,d) = nod%bcv(n,d)
        print *, "nod ", n, ", ",nod%bcv(n,d), ", d", d
      end if
      
      if (nod%is_fix(n,d) .eqv. .true.) then
        nod%v(n,d) = 0.0
      end if 
      d = d + 1 
    end do !dim
    n = n + 1
  end do !Node    
end subroutine

subroutine impose_bca
  implicit none
  integer :: n, d
  n = 1
  do while (n .le. node_count)    
    d = 1
    do while (d .le. 2)
      ! if (nod%is_bcv(n,d) .eqv. .true.) then
        ! nod%v(n,d) = nod%bcv(n,d)
        ! print *, "nod ", n, ", ",nod%bcv(n,d), ", d", d
      ! end if
      
      if (nod%is_fix(n,d) .eqv. .true.) then
        nod%a(n,d) = 0.0
      end if 
      d = d + 1 
    end do !dim
    n = n + 1
  end do !Node    
end subroutine

end module