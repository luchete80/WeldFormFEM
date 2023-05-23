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
  invmat(1,:) = 1.0d0/(det(a))*[ a(2,2),-a(1,2)]
  invmat(2,:) = 1.0d0/(det(a))*[-a(2,1), a(1,1)]
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
  
  integer :: i,j,k
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
    do while (i<=2)
      do while (j<=2)
        !TODO: DO THIS ONCE AT THE BEGINING ONLY ONCE FOR EACH ELEMENT TYPE
        dHrs(1,:)=[(1+s(j)),-(1+s(j)),-(1-s(j)),(1-s(j))]
        dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]   
        dHrs(:,:) = dHrs(:,:)*0.25

        print *, "dHrs ", dHrs
        test = matmul(dHrs,x2)
        print *, "x2, ", x2
        elem%jacob(e,:,:) = test
        elem%dHxy(e,:,:) = matmul(invmat(test),dHrs) !Bathe 5.25
        print *, "inv mat", elem%dHxy(e,:,:)
        
        !DERIVATIVE MATRICES
        !TODO: CHANGE IF DIM != 2
        k=1
        do while (k<nodxelem)
          elem%bl(e,1,dim*(k-1)+k  ) = elem%dHxy(e,1,k)
          elem%bl(e,2,dim*(k-1)+k+1) = elem%dHxy(e,2,k)
          elem%bl(e,3,dim*(k-1)+k  ) = elem%dHxy(e,2,k) 
          elem%bl(e,3,dim*(k-1)+k+1) = elem%dHxy(e,1,k)     

          elem%bnl(e,1,dim*(k-1)+k  ) = elem%dHxy(e,1,k)
          elem%bnl(e,2,dim*(k-1)+k  ) = elem%dHxy(e,2,k)
          elem%bnl(e,3,dim*(k-1)+k+1) = elem%dHxy(e,1,k) 
          elem%bnl(e,4,dim*(k-1)+k+1) = elem%dHxy(e,2,k)     
          k = k+1
        end do
        print *, "jacob e ", elem%jacob(e,:,:)
        
        elem%detJ(e) = det(elem%jacob(e,:,:))
        print *, "det J", elem%detJ(e)
        !print *, "bl ", elem%bl(e,:,:)
        !TODO CHANGE ZERO

        if (dim .eq. 2) then
          temph(1,:) = 0.25*[(1+r(i))*(1+s(j)),0.0d0,(1.0-r(i))*(1+s(j)),0.0d0,(1-r(i))*(1-s(j)),0.0d0,(1+r(i))*(1-s(j)),0.0d0]
          k = 1
          do while (k <= nodxelem)
            temph(2,2*k) = temph(1,2*k-1) !TODO: CHANGE IN 3D
            k = k + 1
          end do
        end if 
        elem%math(e,:,:) = elem%math(e,:,:) + temph(:,:)*elem%detJ(e)
        print *, "mat h ", elem%math(e,:,:)
        !print *, "BL ", elem%bl
        elem%matknl(e,:,:) = elem%matknl(e,:,:) + matmul(matmul(transpose(elem%bnl(e,:,:)),elem%tau(e,:,:)),&
                              &elem%bnl(e,:,:))*elem%detJ(e) !Nodal Weight mat
        elem%matkl(e,:,:) = elem%matkl(e,:,:) + matmul(matmul(transpose(elem%bl(e,:,:)),mat_C),elem%bl(e,:,:))*elem%detJ(e) !Nodal Weight mat
        elem%matm (e,:,:) = elem%matm (e,:,:) + matmul(transpose(elem%math(e,:,:)),elem%math(e,:,:)) *elem%detJ(e)!Mass matrix
        ! print *, "element mat m ", elem%matm (e,:,:)
        j = j +1
      end do
      i = i +1
    end do
    
    print *, "element mat m ", elem%matm (e,:,:)
    
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

!NEEDED FOR STRAIN AND INTERNAL FORCES CALCULATION
!IS REALLY NEEDED TO STORE?
subroutine disassemble_uele()
  integer :: e, i, n
  do while (e .le. elem_count)
    !print *, "elem ", e
    n = 1
    do while (n .le. nodxelem)
      do while (i .le. dim )
        elem%uele (e, 2*n-1,1) = uglob(2*elem%elnod(e,n)-1,1)
        i = i + 1
      end do
      n = n + 1
    end do ! Element node
    e = e + 1
  end do ! e
end subroutine

subroutine assemble_int_forces()
  integer :: e, i, n, iglob
  real(fp_kind), dimension(nodxelem*dim,1) :: utemp, rtemp
  
  print *, "assemblying int forces"
  rint_glob (:,:) = 0.0d0
  e = 1
  do while (e .le. elem_count)
    !print *, "elem ", e
    rtemp = matmul(elem%matkl(e,:,:) + elem%matknl(e,:,:), elem%uele(e,:,:))
    n = 1
    do while (n .le. nodxelem)
      i = 1
      print *,"elem mat kl", elem%matkl(e,:,:)
      do while (i .le. dim )
        iglob  = dim * (elem%elnod(e,n) - 1 ) + i
        rint_glob(elem%elnod(e,n),i) =  rint_glob(elem%elnod(e,n),i) + rtemp(dim*(n-1)+i,1)
        i = i + 1
      end do !element row
      n = n + 1
    end do ! Element node
    e = e + 1
  end do ! e
end subroutine 

end module Matrices