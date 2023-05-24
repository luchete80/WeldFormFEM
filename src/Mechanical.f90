module Mechanical
use Domain
use ElementData

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

end module