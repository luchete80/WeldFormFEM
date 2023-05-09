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
  real(fp_kind), dimension(dim, dim*nodxelem) :: temph
  
  integer :: i,j
  real(fp_kind), dimension(2) :: r, s
  e = 1
  r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
  s(1) = r(1)          ; s(2) =  r(2)
  !! Update x2 vector (this is useful for strain and stress things)
  
  
  do while (e <= elem_count)
    print *, "el ", e 
    
    elem%matkl(e,:,:) = 0.0
    elem%matknl(e,:,:) = 0.0
    elem%matm(e,:,:) = 0.0
    
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
          elem%bl(e,1,dim*(i-1)+i  ) = elem%dHxy(e,1,i)
          elem%bl(e,2,dim*(i-1)+i+1) = elem%dHxy(e,2,i)
          elem%bl(e,3,dim*(i-1)+i  )   = elem%dHxy(e,2,i) 
          elem%bl(e,3,dim*(i-1)+i+1) = elem%dHxy(e,1,i)     

          elem%bnl(e,1,dim*(i-1)+i  ) = elem%dHxy(e,1,i)
          elem%bnl(e,2,dim*(i-1)+i  ) = elem%dHxy(e,2,i)
          elem%bnl(e,3,dim*(i-1)+i+1) = elem%dHxy(e,1,i) 
          elem%bnl(e,4,dim*(i-1)+i+1) = elem%dHxy(e,2,i)     
          i= i+1
        end do
        print *, "jacob e ", elem%jacob(e,:,:)
        print *, "det J", elem%detJ(e)
        elem%detJ(e) = det(elem%jacob(e,:,:))
        !TODO CHANGE ZERO
        if (dim .eq. 2) then
          temph(1,:) = 0.2*[(1+r(i))*(1+s(j)),0.0d0,(1.0-r(i))*(1+s(j)),0.0d0,(1-r(i))*(1-s(j)),0.0d0,(1+r(i))*(1-s(j)),0.0d0]
          i = 1
          do while (i < nodxelem)
            temph(2,2*i) = temph(1,2*i-1)
            i = i + 1
          end do
        end if 
        elem%math(e,:,:) = elem%math(e,:,:) + temph(:,:)*elem%detJ(e)
        print *, "element mat m ", elem%math (e,:,:)
        !print *, "BL ", elem%bl
        elem%matkl(e,:,:) = elem%matkl(e,:,:) + matmul(matmul(transpose(elem%bl(e,:,:)),mat_C),elem%bl(e,:,:))*elem%detJ(e) !Nodal Weight mat
        elem%matm (e,:,:) = elem%matm (e,:,:) + matmul(transpose(elem%math(e,:,:)),elem%math(e,:,:)) !Mass matrix
        print *, "element mat m ", elem%matm (e,:,:)
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

subroutine assemble_mass_matrix ()
  integer :: e, i, j, n, iglob, jglob
  
  m_glob (:,:) = 0.0d0
  e = 1
  do while (e .le. elem_count)
	n = 1
	do while (n .le. nodxelem)
		i = 1
		do while (i .le. dim )
			j = 1
			do while (j .le. dim )
				iglob  = dim * (elem%elnod(e,n) - 1 ) + i
				jglob  = dim * (elem%elnod(e,n) - 1 ) + j
				m_glob(iglob,jglob) = m_glob(iglob,jglob) + elem%matm (e,i,j)
				j = j + 1
			end do
			i = i + 1
		end do !element row
		n = n + 1
	end do ! Element node
    e = e + 1
  end do ! e
end subroutine

subroutine assemble_int_forces()
  integer :: e, i, j, n, iglob
  real(fp_kind), dimension(nodxelem*dim,1) :: utemp, rtemp
  
  print *, "assemblying int forces"
  rint_glob (:,:) = 0.0d0
  e = 1
  do while (e .le. elem_count)
    !print *, "elem ", e
  n = 1
	do while (n .le. nodxelem)
		i = 1
		do while (i .le. dim )
      iglob  = dim * (elem%elnod(e,n) - 1 ) + i
      rtemp = matmul(elem%matkl(e,:,:), utemp(:,:))
      rint_glob(elem%elnod(e,n),i) =  rint_glob(elem%elnod(e,n),i) + rtemp(dim*(n-1)+i,1)
			i = i + 1
		end do !element row
		n = n + 1
	end do ! Element node
    e = e + 1
  end do ! e
end subroutine 

end module Matrices