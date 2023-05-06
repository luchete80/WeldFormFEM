module Matrices
use Domain
use ElementData

implicit none 

contains 
subroutine calculate_element_matrices ()
  integer :: e
  ! !rg=gauss[ig]
  ! !sg=gauss[jg]
  real(fp_kind), dimension(dim,nodxelem) :: dHrs
  real(fp_kind), dimension(nodxelem,dim) :: x2
  real(fp_kind), dimension(dim,dim) :: test
  
  integer :: i,j
  real(fp_kind) :: r, s
  e = 1
  
  !! Update x2 vector (this is useful for strain and stress things)
  
  do while (e < elem_count)
    i=1
    ! do while (i.le.nodxelem)
        ! x2(:,:)=nod%x(i, elem%elnod(e,i))
    ! end do
    ! This could be done once
    do while (i<2)
      do while (j<2)
        !TODO: DO THIS ONCE AT THE BEGINING
        dHrs(1,:)=[-(1-s),(1-s),-(1+s),(1+s)]
        dHrs(2,:)=[-(1-r),-(1+r),(1-r),(1+r)]   
        dHrs(:,:) = dHrs(:,:)*0.25
      end do
    end do
    test = matmul(dHrs,x2)
    elem%jacob(e,:,:) = test
    !elem%dHxy(e,:,:) = matmul(inv(test),dHrs)
    
    !For node
    i=1
    do while (i<nodxelem)
      elem%bl(e,1,2*i  ) = elem%dHxy(e,1,i)
      elem%bl(e,2,2*i+1) = elem%dHxy(e,2,i)
      elem%bl(e,3,2*i)   = elem%dHxy(e,2,i) 
      elem%bl(e,3,2*i+1) = elem%dHxy(e,1,i)     

      elem%bnl(e,1,2*i  ) = elem%dHxy(e,1,i)
      elem%bnl(e,2,2*i  ) = elem%dHxy(e,2,i)
      elem%bnl(e,3,2*i+1) = elem%dHxy(e,1,i) 
      elem%bnl(e,4,2*i+1) = elem%dHxy(e,2,i)         
    end do
    e = e + 1 
    ! #Numerated as in Bathe
    ! Ns  =0.25*matrix([(1+sg)*(1+rg),(1-rg)*(1+sg),(1-sg)*(1-rg),(1-sg)*(1+rg)])   
    ! dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
    ! #Numerated as in deal.ii
    ! #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
    ! dHrs/=4.
    ! J=dHrs*X2
    ! dHxy=linalg.inv(J)*dHrs
    ! detJ=linalg.det(J)
  end do
end subroutine

end module Matrices