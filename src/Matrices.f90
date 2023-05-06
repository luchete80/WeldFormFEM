module Matrices
use Domain
use ElementData

implicit none 

contains 

function det (a)
  real(fp_kind), dimension(dim,dim), intent (in) :: a 
  real(fp_kind) :: det
  !if (dim .eq. 2) then
  det = a(1,1)*a(2,2)-a(1,2)*a(2,1)
  !end if
end function

function invmat (a)
  real(fp_kind), dimension(dim,dim), intent (in) :: a 
  real(fp_kind), dimension(dim,dim) :: invmat 
  !if (dim .eq. 2) then
  invmat(1,:) = [ a(2,2),-a(1,2)]
  invmat(2,:) = [-a(2,1), a(1,1)]
  !end if
end function

subroutine calculate_element_matrices ()
  integer :: e
  ! !rg=gauss[ig]
  ! !sg=gauss[jg]
  real(fp_kind), dimension(dim,nodxelem) :: dHrs
  real(fp_kind), dimension(nodxelem,dim) :: x2
  real(fp_kind), dimension(dim,dim) :: test
  
  integer :: i,j
  real(fp_kind), dimension(2) :: r, s
  e = 1
  r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
  s(1) = r(1)          ; s(2) =  r(2)
  !! Update x2 vector (this is useful for strain and stress things)
  
  do while (e < elem_count)
    print *, "el ", e
    i=1
    print *, "x2 ", x2
    do while (i.le.nodxelem)
        x2(i,:)=nod%x(elem%elnod(e,i),:)
    end do
    ! TODO: This could be done once
    i = 1; j = 1
    do while (i<2)
      do while (j<2)
        !TODO: DO THIS ONCE AT THE BEGINING
        dHrs(1,:)=[-(1-s(j)),(1-s(j)),-(1+s(j)),(1+s(j))]
        dHrs(2,:)=[-(1-r(i)),-(1+r(i)),(1-r(i)),(1+r(i))]   
        dHrs(:,:) = dHrs(:,:)*0.25
      end do
    end do
    test = matmul(dHrs,x2)
    elem%jacob(e,:,:) = test
    elem%dHxy(e,:,:) = matmul(invmat(test),dHrs)
    
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
    elem%detJ(e) = det(elem%jacob(e,:,:))
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