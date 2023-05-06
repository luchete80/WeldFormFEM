module Matrices
use Domain
use ElementData

implicit none 

contains 

! subroutine printarray(Grid)
  ! implicit none
  ! CHARACTER(len=1) :: Grid(3,2)
  ! integer :: i
  ! Grid = "x"
  ! do i = 1, ubound(Grid, 1)
     ! print *, Grid(i, :)
  ! end do
! end subroutine

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
  real(fp_kind), dimension(dim*nodxelem,dim*nodxelem) ::tempk
  
  integer :: i,j
  real(fp_kind), dimension(2) :: r, s
  e = 1
  r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
  s(1) = r(1)          ; s(2) =  r(2)
  !! Update x2 vector (this is useful for strain and stress things)
  
  
  do while (e < elem_count)
    print *, "el ", e
    
    elem%matkl(e,:,:) = 0.0
    elem%matknl(e,:,:) = 0.0
    
    print *, "nodxelem ", nodxelem
    i=1
    do while (i.le.nodxelem)
        print *, "elnod " , elem%elnod(e,i)
        x2(i,:)=nod%x(elem%elnod(e,i),:)
        i = i+1
    end do
    !print *, "x2 ", x2
    !printarray(x2)
    ! TODO: This could be done once
    i = 1; j = 1 !TODO: CHANGE TO PLAIN DO (IN ORDER TO INCLUDE 3D)
    do while (i<2)
      do while (j<2)
        !TODO: DO THIS ONCE AT THE BEGINING
        dHrs(1,:)=[-(1-s(j)),(1-s(j)),-(1+s(j)),(1+s(j))]
        dHrs(2,:)=[-(1-r(i)),-(1+r(i)),(1-r(i)),(1+r(i))]   
        dHrs(:,:) = dHrs(:,:)*0.25

        print *, "dHrs ", dHrs
        test = matmul(dHrs,x2)
        elem%jacob(e,:,:) = test
        elem%dHxy(e,:,:) = matmul(invmat(test),dHrs)
        
        !DERIVATIVE MATRICES
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
          i= i+1
        end do
        elem%detJ(e) = det(elem%jacob(e,:,:))
        
        print *, "BL ", elem%bl
        elem%matkl(e,:,:) = elem%matkl(e,:,:) + matmul(matmul(elem%bl(e,:,:),mat_C),elem%bl(e,:,:))*elem%detJ(e)
        j = j +1
      end do
      i = i +1
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