module Mechanical
use Domain
use ElementData

subroutine calculate_element_strains ()
  integer :: e
  e = 1
  
  integer :: i,j,k
  real(fp_kind), dimension(2) :: r, s
  
  r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
  s(1) = r(1)          ; s(2) =  r(2)
  do while (e <= elem_count)  
    !Is only linear matrix?
    
    elem%eps(e,:,:) = elem%bl(e,

    e = e + 1
  end do
end subroutine