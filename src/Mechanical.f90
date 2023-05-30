module Mechanical
use Domain
use ElementData
use NodeData

contains
!THIS SHOULD BE DONE AT t+1/2dt
subroutine cal_elem_strains ()
  implicit none
  integer :: e, i,j,k, gp
  real(fp_kind), dimension(2) :: r, s
  
 
  r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
  s(1) = r(1)          ; s(2) =  r(2)

  gp = 1
  do e=1, elem_count
    !Is only linear matrix?    
    elem%strain(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%uele (e,:,:)) 

  end do
end subroutine

subroutine calc_elem_vol ()
  implicit none
  integer :: e
  real(fp_kind), dimension(2) :: r, s
  r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
  s(1) = r(1)          ; s(2) =  r(2)
  
  ! P00+(Cs0*Cs0)*(Density-Density0);
  do e = 1, elem_count
    elem%vol(e) = 0.0d0
    !elem%vol(e) = 
    !elem%vol(e) = elem%vol(e) + 
  
  end do

end subroutine

!!!!!EOS: Equation of State
subroutine calc_elem_density ()
  implicit none
  integer :: e
  ! P00+(Cs0*Cs0)*(Density-Density0);
  do e = 1, elem_count
    !elem%rho(e) = 
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