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

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

subroutine M33INV (A, AINV, OK_FLAG)

  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
  LOGICAL, INTENT(OUT) :: OK_FLAG

  DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
  DOUBLE PRECISION :: DETE
  DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


  ! DET =   A(1,1)*A(2,2)*A(3,3)  &
        ! - A(1,1)*A(2,3)*A(3,2)  &
        ! - A(1,2)*A(2,1)*A(3,3)  &
        ! + A(1,2)*A(2,3)*A(3,1)  &
        ! + A(1,3)*A(2,1)*A(3,2)  &
        ! - A(1,3)*A(2,2)*A(3,1)
  
  dete = det(A)

  IF (ABS(DETE) .LE. EPS) THEN
     AINV = 0.0D0
     OK_FLAG = .FALSE.
     RETURN
  END IF

  COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
  COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
  COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
  COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
  COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

  AINV = TRANSPOSE(COFACTOR) / DETE

  OK_FLAG = .TRUE.

RETURN

end subroutine M33INV

function det (a)
  real(fp_kind), dimension(dim,dim), intent (in) :: a 
  real(fp_kind) :: det
  if (dim .eq. 2) then
    det = a(1,1)*a(2,2)-a(1,2)*a(2,1)
  else 
  DET =   A(1,1)*A(2,2)*A(3,3)  &
        - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3)  &
        + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2)  &
        - A(1,3)*A(2,2)*A(3,1)  
  end if
end function

function invmat (a)
  real(fp_kind), dimension(dim,dim), intent (in) :: a 
  real(fp_kind), dimension(dim,dim) :: invmat 
  !if (dim .eq. 2) then
  invmat(1,:) = 1.0d0/(det(a))*[ a(2,2),-a(1,2)]
  invmat(2,:) = 1.0d0/(det(a))*[-a(2,1), a(1,1)]
  !end if
end function

function adj (a)
  real(fp_kind), dimension(dim,dim), intent (in) :: a 
  real(fp_kind), dimension(dim,dim) :: cofactor,adj
  
  if (dim .eq. 2) then
    adj(1,:) = [ a(2,2),-a(1,2)]
    adj(2,:) = [-a(2,1), a(1,1)]
  else
    cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    
    adj = TRANSPOSE(COFACTOR)
  end if
end function


!!!! IN ORDER TO CALCULATE IT ONLY ONCE
!!!!!!!!!!!! IN RS PLANE ---->>>>
!!!!!!!! 4------3
!!!!!!!! |      |
!!!!!!!! 1 ---- 2
!!!! CALCULATE JACOBIAN AND DETERMINANT
subroutine calculate_element_Jacobian ()
  integer :: e
  ! !rg=gauss[ig]
  ! !sg=gauss[jg]
  real(fp_kind), dimension(dim,nodxelem) :: dHrs !!! USED ONLY FOR SEVERAL GAUSS POINTS
  real(fp_kind), dimension(nodxelem,dim) :: x2
  real(fp_kind), dimension(dim,dim) :: test
  real(fp_kind), dimension(dim, dim*nodxelem) :: temph
  
  integer :: i,j,k, gp
  real(fp_kind):: r   !!! USED ONLY FOR SEVERAL GAUSS POINTS
  real(fp_kind), dimension(8,3):: gpc !!! gauss point coordinates, r,s,t
  
  gp = 1
  do e=1, elem_count
    print *, "el ", e 
    
    do i=1,nodxelem
        !print *, "elnod " , elem%elnod(e,i)
        x2(i,:)=nod%x(elem%elnod(e,i),:)
    end do
    
    if (elem%gausspc(e) .eq. 1) then      
    
      if (dim .eq. 2) then 
      
        else !!!DIM 3
          !!!!! SETTING LIKE THIS AVOID MATMUL
          elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)-x2(5,:)+x2(6,:)+x2(7,:)-x2(8,:)
          elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)-x2(5,:)-x2(6,:)+x2(7,:)+x2(8,:)
          elem%jacob(e,gp,3,:) = -x2(1,:)-x2(2,:)-x2(3,:)-x2(4,:)+x2(5,:)+x2(6,:)+x2(7,:)+x2(8,:)
          !elem%jacob(e,gp,2,:) = [-x2(1,2),-x2(2,2), x2(3,2), x2(4,2),-x2(5,2),-x2(6,2), x2(7,2), x2(8,2)]
          !elem%jacob(e,gp,3,:) = [-x2(1,3),-x2(2,3), x2(3,3), x2(4,3),-x2(5,3),-x2(6,3), x2(7,3), x2(8,3)]
          ! dHrs(1,:)=[-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0] AND THIS IS dHrs*x2
          ! dHrs(2,:)=[-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0]       
          ! dHrs(3,:)=[-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0]  
          ! elem%jacob(e,gp,1,:) = matmul(dHrs,x2)
      end if  !!!!DIM
      elem%jacob(e,gp,:,:) = 0.125*elem%jacob(e,gp,:,:)
    else !!!!! GP > 1
      r = 1.0/sqrt(3.0);
      gpc(1,:)=[-r,-r,-r];   gpc(2,:)=[ r,-r,-r];      gpc(3,:)=[-r, r,-r];      gpc(4,:)=[ r, r,-r];
      gpc(1,:)=[-r,-r, r];   gpc(2,:)=[ r,-r, r];      gpc(3,:)=[-r, r, r];      gpc(4,:)=[ r, r, r];
      do i=1,nodxelem
          x2(i,:)=nod%x(elem%elnod(e,i),:)
      end do
    
      if (dim .eq. 3) then
        do gp = 1,8
          dHrs(1,:)=[-1.0*(1-gpc(gp,2))*(1.0-gpc(gp,3)),(1-gpc(gp,2))*(1.0-gpc(gp,3))&
                          ,(1+gpc(gp,2))*(1.0-gpc(gp,2)),-1.0*(1+gpc(gp,2))*(1.0-gpc(gp,2)) &
                    ,-1.0*(1-gpc(gp,2))*(1.0-gpc(gp,3)),(1-gpc(gp,2))*(1.0-gpc(gp,3))&
                    ,(1+gpc(gp,2))*(1.0-gpc(gp,2)),-1.0*(1+gpc(gp,2))*(1.0-gpc(gp,2))]
          !dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]         
          !dHrs(3,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))] 
          elem%jacob(e,gp,:,:) = matmul(dHrs,x2)
        end do !gp
      end if
    end if !!gp ==1
    elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
  end do !element

end subroutine

!!!!! PREVIOUSLY JACOBIAN DETERMINANT SHOULD BE CALCULATED
subroutine calculate_element_shapeMat ()
  integer :: e,d
  ! !rg=gauss[ig]
  ! !sg=gauss[jg]
  real(fp_kind), dimension(dim,nodxelem) :: dHrs
  real(fp_kind), dimension(nodxelem,dim) :: x2
  real(fp_kind), dimension(dim,dim) :: test
  real(fp_kind), dimension(dim, dim*nodxelem) :: temph
  
  integer :: i,j,k, gp
  real(fp_kind):: r
  real(fp_kind), dimension(8,3):: gpc !!! gauss point coordinates, r,s,t
  !! Update x2 vector (this is useful for strain and stress things)
  

  do e=1, elem_count
    gp = 1
    if (elem%gausspc(e) .eq. 1) then
      elem%math(e,gp,:,:) = 0.0d0
      if (dim .eq. 2) then 
      
        else !!!DIM 3
          temph(:,:) = 0.0d0
          do k =1, nodxelem
              !temph(d,dim*(k-1)+d) = 0.125 !TODO: CHANGE IN 3D
              elem%math(e,gp,:,k)=0.125
              !print *, "temp h i ", d, "j ", dim*(k-1)+d, ": "
          end do
          ! print *, "temp h", temph(:,:)
          ! temph(1,:) = [1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0]
          ! temph(2,:) = [0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0]
          ! temph(3,:) = [0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0]

      end if  !dim
        
        elem%matm(e,:,:) = matmul(transpose(elem%math(e,gp,:,:)),elem%math(e,gp,:,:))*elem%rho(e)*elem%detJ(e,gp)*8.0 !!!2.0 ^3 WEIGHT
        
        !print *, "MAT M", elem%matm(e,:,:)
      
    else !!!!! GP > 1
      if (dim .eq. 2) then 
      
        else !!!DIM 3
          gpc(1,:)=[-r,-r,-r];   gpc(2,:)=[ r,-r,-r];      gpc(3,:)=[-r, r,-r];      gpc(4,:)=[ r, r,-r];
          gpc(5,:)=[-r,-r, r];   gpc(6,:)=[ r,-r, r];      gpc(7,:)=[-r, r, r];      gpc(8,:)=[ r, r, r];
          do gp = 1,8
            elem%math(e,gp, 1,:) = 0.125*[(1-gpc(gp,1))*(1-gpc(gp,2))*(1-gpc(gp,3)),(1+gpc(gp,1))*(1-gpc(gp,2))*(1-gpc(gp,3)), &
                                (1+gpc(gp,1))*(1+gpc(gp,2))*(1-gpc(gp,3)),(1-gpc(gp,1))*(1+gpc(gp,2))*(1-gpc(gp,3)), &
                                (1-gpc(gp,1))*(1-gpc(gp,2))*(1+gpc(gp,3)),(1+gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3)), &
                                (1+gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3)),(1-gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3))]
            !elem%matm(e,:,:) = elem%matm(e,:,:) + matmul(transpose(elem%math(e,gp,:,:)),elem%math(e,gp,:,:))*elem%rho(e)*elem%detJ(e,gp)*8.0 !!!2.0 ^3 WEIGHT
          end do
            ! k = 1
            ! do while (k <= nodxelem)
              ! temph(2,2*k) = temph(1,2*k-1) !TODO: CHANGE IN 3D
              ! k = k + 1
            ! end do
 
          ! elem%math(e,gp,:,:) = temph(:,:)*elem%detJ(e,gp)
          !print *, "mat h ", elem%math(e,gp,:,:)        
      end if !!! if dim
    end if ! GP>1
  end do !element
end subroutine

!!!!!dxdxy and B matrices
!!!!! WITHOUT CALCULATING DETERMINANT, ONLY ADJOINT
subroutine calculate_element_derivMat ()
  integer :: e,d
  ! !rg=gauss[ig]
  ! !sg=gauss[jg]
  real(fp_kind), dimension(dim,nodxelem) :: dHrs
  real(fp_kind), dimension(dim,dim) :: test, invJ
  
  integer :: i,j,k, gp
  real(fp_kind), dimension(2) :: r, s

  !! Update x2 vector (this is useful for strain and stress things)
  

  do e=1, elem_count
    !print *, "el ", e 

    gp = 1
    if (elem%gausspc(e) .eq. 1) then
    
      if (dim .eq. 2) then       
        else !!!DIM 3
          ! dHrs(1,:)=[-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0]
          ! dHrs(2,:)=[-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0]       
          ! dHrs(3,:)=[-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0]  
          ! elem%jacob(e,gp,:,:) = matmul(dHrs,x2)
          !elem%dHxy(e,gp,:,:) = matmul(invmat(test),dHrs) !Bathe 5.25
          invJ = adj(elem%jacob(e,gp,:,:))/elem%detJ(e,gp) !!!! ALREADY CALCULATED         
          !elem%dHxy(e,gp,:,:) = matmul(invmat(test),dHrs) !Bathe 5.25
          !!!! DONE LIKE THIS TO AVOID MULTS
          
          !!!! REMAINS ELEM_VOLUME
          elem%dHxy(e,gp,:,1) = -invJ(:,1)-invJ(:,2)-invJ(:,3) !For each 3 rows of inv J and dHdxy
          elem%dHxy(e,gp,:,2) =  invJ(:,1)-invJ(:,2)-invJ(:,3)
          elem%dHxy(e,gp,:,3) =  invJ(:,1)+invJ(:,2)-invJ(:,3)
          elem%dHxy(e,gp,:,4) = -invJ(:,1)+invJ(:,2)-invJ(:,3)
          elem%dHxy(e,gp,:,5) = -invJ(:,1)-invJ(:,2)+invJ(:,3)
          elem%dHxy(e,gp,:,6) =  invJ(:,1)-invJ(:,2)+invJ(:,3)
          elem%dHxy(e,gp,:,7) =  invJ(:,1)+invJ(:,2)+invJ(:,3)
          elem%dHxy(e,gp,:,8) = -invJ(:,1)+invJ(:,2)+invJ(:,3)
          
          ! print *, "1,1", elem%dHxy(e,gp,1,1), "inv 1,1 1,2 1,3", invJ(1,1), invJ(1,2), invJ(1,3)
          ! print *, "1,1", elem%dHxy(e,gp,2,1), "inv 1,1 1,2 1,3", invJ(2,1), invJ(2,2), invJ(2,3)
          do k=1,nodxelem
            do d=1,dim
            elem%bl(e,gp,d,dim*(k-1)+d ) = elem%dHxy(e,gp,d,k) 
            end do
            elem%bl(e,gp,4,dim*(k-1)+1) = elem%dHxy(e,gp,2,k)
            elem%bl(e,gp,4,dim*(k-1)+2) = elem%dHxy(e,gp,1,k)
            
            elem%bl(e,gp,5,dim*(k-1)+2) = elem%dHxy(e,gp,3,k)
            elem%bl(e,gp,5,dim*(k-1)+3) = elem%dHxy(e,gp,2,k)
            
            elem%bl(e,gp,6,dim*(k-1)+1) = elem%dHxy(e,gp,3,k)
            elem%bl(e,gp,6,dim*(k-1)+3) = elem%dHxy(e,gp,1,k)         
          end do
      end if  !!!!DIM
      
    else !!!!! GP > 1

      r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
      s(1) = r(1)          ; s(2) =  r(2)      
      i = 1; j = 1 !TODO: CHANGE TO PLAIN DO (IN ORDER TO INCLUDE 3D)
      do while (i<=2)
        j = 1
        do while (j<=2)
          !TODO: DO THIS ONCE AT THE BEGINING ONLY ONCE FOR EACH ELEMENT TYPE
          dHrs(1,:)=[(1+s(j)),-(1+s(j)),-(1-s(j)),(1-s(j))]
          dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]   
          dHrs(:,:) = dHrs(:,:)*0.25
          elem%dHxy(e,gp,:,:) = matmul(invmat(test),dHrs) !Bathe 5.25
          !print *, "inv mat", elem%dHxy(e,gp,:,:)
          
          !DERIVATIVE MATRICES
          !TODO: CHANGE IF DIM != 2
          k=1
          do while (k<nodxelem)
            elem%bl(e,gp,1,dim*(k-1)+k  ) = elem%dHxy(e,gp,1,k)
            elem%bl(e,gp,2,dim*(k-1)+k+1) = elem%dHxy(e,gp,2,k)
            elem%bl(e,gp,3,dim*(k-1)+k  ) = elem%dHxy(e,gp,2,k) 
            elem%bl(e,gp,3,dim*(k-1)+k+1) = elem%dHxy(e,gp,1,k)     
            k = k+1
          end do
          ! print *, "element mat m ", elem%matm (e,:,:)
          gp = gp + 1
          j = j +1
        end do
        i = i +1
        
        elem%matm (e,:,:) = elem%matm (e,:,:) * elem%rho(e) !!ASUMMING constant element density
      end do !i
    end if   !j
    !print *, "element",e," mat m ", elem%matm (e,:,:)

  end do
end subroutine

!!!!!!!!!!!!!!!!!!!!
!!!!! OLD 
!!!! ---------
subroutine calculate_element_matrices ()
  integer :: e
  ! !rg=gauss[ig]
  ! !sg=gauss[jg]
  real(fp_kind), dimension(dim,nodxelem) :: dHrs
  real(fp_kind), dimension(nodxelem,dim) :: x2
  real(fp_kind), dimension(dim,dim) :: test
  real(fp_kind), dimension(dim, dim*nodxelem) :: temph
  
  integer :: i,j,k, gp
  real(fp_kind), dimension(2) :: r, s

  !! Update x2 vector (this is useful for strain and stress things)
  

  do e=1, elem_count
    print *, "el ", e 
    
    elem%matkl(e,:,:) = 0.0
    elem%matknl(e,:,:) = 0.0
    if (calc_m .eqv. .True.) then
      elem%matm(e,:,:) = 0.0
    end if 
    
    print *, "nodxelem ", nodxelem
    
    do i=1,nodxelem
        !print *, "elnod " , elem%elnod(e,i)
        x2(i,:)=nod%x(elem%elnod(e,i),:)
    end do
    !print *, "x2 ", x2
    !printarray(x2)
    ! TODO: This could be done once
    gp = 1
    if (elem%gausspc(e) .eq. 1) then
    
      if (dim .eq. 2) then 
      
        else !!!DIM 3
          dHrs(1,:)=[-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0]
          dHrs(2,:)=[-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0]       
          dHrs(3,:)=[-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0]  
          elem%jacob(e,gp,:,:) = matmul(dHrs,x2)
      end if  !!!!DIM
      
    else !!!!! GP > 1

      r(1) = -1.0/sqrt(3.0); r(2) = -r(1)
      s(1) = r(1)          ; s(2) =  r(2)      
      i = 1; j = 1 !TODO: CHANGE TO PLAIN DO (IN ORDER TO INCLUDE 3D)
      do while (i<=2)
        j = 1
        do while (j<=2)
          !TODO: DO THIS ONCE AT THE BEGINING ONLY ONCE FOR EACH ELEMENT TYPE
          dHrs(1,:)=[(1+s(j)),-(1+s(j)),-(1-s(j)),(1-s(j))]
          dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]   
          dHrs(:,:) = dHrs(:,:)*0.25

          !print *, "dHrs ", dHrs
          test = matmul(dHrs,x2)
          !print *, "x2, ", x2
          elem%jacob(e,gp,:,:) = test
          elem%dHxy(e,gp,:,:) = matmul(invmat(test),dHrs) !Bathe 5.25
          !print *, "inv mat", elem%dHxy(e,gp,:,:)
          
          !DERIVATIVE MATRICES
          !TODO: CHANGE IF DIM != 2
          k=1
          do while (k<nodxelem)
            elem%bl(e,gp,1,dim*(k-1)+k  ) = elem%dHxy(e,gp,1,k)
            elem%bl(e,gp,2,dim*(k-1)+k+1) = elem%dHxy(e,gp,2,k)
            elem%bl(e,gp,3,dim*(k-1)+k  ) = elem%dHxy(e,gp,2,k) 
            elem%bl(e,gp,3,dim*(k-1)+k+1) = elem%dHxy(e,gp,1,k)     
            k = k+1
          end do
          print *, "jacob e ", elem%jacob(e,gp,:,:)
          
          elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
          !print *, "det J", elem%detJ(e,gp)
          !print *, "bl ", elem%bl(e,gp,:,:)
          !TODO CHANGE ZERO

          if (dim .eq. 2) then
            temph(1,:) = 0.25*[(1+r(i))*(1+s(j)),0.0d0,(1.0-r(i))*(1+s(j)),0.0d0,(1-r(i))*(1-s(j)),0.0d0,(1+r(i))*(1-s(j)),0.0d0]
            k = 1
            do while (k <= nodxelem)
              temph(2,2*k) = temph(1,2*k-1) !TODO: CHANGE IN 3D
              k = k + 1
            end do
          end if 
          elem%math(e,gp,:,:) = temph(:,:)*elem%detJ(e,gp)
          !print *, "mat h ", elem%math(e,gp,:,:)
          !print *, "BL ", elem%bl
          ! elem%matknl(e,:,:) = elem%matknl(e,:,:) + matmul(matmul(transpose(elem%bnl(e,gp,:,:)),elem%tau(e,gp,:,:)),&
                                ! &elem%bnl(e,gp,:,:))*elem%detJ(e,gp) !Nodal Weight mat
          ! elem%matkl(e,:,:) = elem%matkl(e,:,:) + matmul(matmul(transpose(elem%bl(e,gp,:,:)),mat_C),elem%bl(e,gp,:,:))*elem%detJ(e,gp) !Nodal Weight mat
          if (calc_m .eqv. .True.) then
            elem%matm (e,:,:) = elem%matm (e,:,:) + matmul(transpose(elem%math(e,gp,:,:)),elem%math(e,gp,:,:)) *elem%detJ(e,gp)!Mass matrix
          end if
          ! print *, "element mat m ", elem%matm (e,:,:)
          gp = gp + 1
          j = j +1
        end do
        i = i +1
        
        elem%matm (e,:,:) = elem%matm (e,:,:) * elem%rho(e) !!ASUMMING constant element density
      end do !i
    end if   !j
    !print *, "element",e," mat m ", elem%matm (e,:,:)

  end do
end subroutine

!!!!! IF MASS MATRIX IS OBTAINED FOR A ZERO FILLED SHAPEMATRIX (dim x dim*nodxelem)
! subroutine assemble_mass_matrix ()
  ! integer :: e,gp, i, j, dof1,dof2, n2
  
  ! m_glob (:,:) = 0.0d0
  ! do e = 1, elem_count
    ! !print *, "elem ", e
    ! do dof1 =1, (nodxelem*dim)
      ! do dof2=1, (nodxelem*dim)
            ! !print *, "elem ", e, "node ", n, " i j matm ",i, j, elem%matm (e,dim*(n-1)+i,dim*(n2-1)+j)            
            ! ! iglob  = dim * (elem%elnod(e,n) - 1 ) + i
            ! ! jglob  = dim * (elem%elnod(e,n2) - 1 ) + j
            ! ! print *, "iloc, jloc ",dim*(n-1)+i, dim*(n2-1)+j, "iglob, jglob", iglob,jglob
            ! m_glob(elem%dof(e,dof1),elem%dof(e,dof2)) = m_glob(elem%dof(e,dof1),elem%dof(e,dof2)) + elem%matm (e,dof1,dof2)
      ! end do !dof1
    ! end do ! dof2 
  ! end do ! e
! end subroutine

subroutine assemble_mass_matrix ()
  integer :: e,gp, i, j, iglob,jglob, n1,n2
  
  m_glob (:,:) = 0.0d0
  do e = 1, elem_count
    !print *, "elem ", e
    do n1 =1, nodxelem
      do n2=1, nodxelem
            !print *, "elem ", e, "node ", n, " i j matm ",i, j, elem%matm (e,dim*(n-1)+i,dim*(n2-1)+j)            
            iglob  = dim * (elem%elnod(e,n1) - 1 ) + i
            jglob  = dim * (elem%elnod(e,n2) - 1 ) + j
            ! print *, "iloc, jloc ",dim*(n-1)+i, dim*(n2-1)+j, "iglob, jglob", iglob,jglob
            m_glob(iglob,jglob) = m_glob(iglob,jglob) + elem%matm (e,n1,n2)
      end do !dof1
    end do ! dof2 
  end do ! e
end subroutine

!NEEDED FOR STRAIN AND INTERNAL FORCES CALCULATION
!IS REALLY NEEDED TO STORE?
!TODO: CHANGE TO DOF
subroutine disassemble_uvele()
  integer :: e, i, n
  do e=1,elem_count
    !print *, "elem ", e
    do n =1,nodxelem
      do i =1, dim 
        !print *, "e ", e, "n ", n, "uele estirado", 2*(n-1)+i , ",global ",elem%elnod(e,n)      
        elem%uele (e,dim*(n-1)+i,1) = nod%u(elem%elnod(e,n),i)
        elem%vele (e,dim*(n-1)+i,1) = nod%v(elem%elnod(e,n),i)
      end do
    end do ! Element node
  end do ! e
end subroutine

subroutine assemble_forces()
  integer :: e, i, n, iglob
  real(fp_kind), dimension(nodxelem*dim,1) :: utemp, rtemp
  
  print *, "assemblying int forces"
  rint_glob (:,:) = 0.0d0
  do e = 1, elem_count
    do n = 1, nodxelem
      !print *,"elem mat kl", elem%matkl(e,:,:)
      !print *, "elem fext ", elem%f_ext(e,n,:)
      do i=1,dim 
        iglob  = dim * (elem%elnod(e,n) - 1 ) + i
        !rint_glob(elem%elnod(e,n),i) =  rint_glob(elem%elnod(e,n),i) + rtemp(dim*(n-1)+i,1)
        rint_glob(elem%elnod(e,n),i) =  rint_glob(elem%elnod(e,n),i) + elem%f_int(e,n,i)
        fext_glob(elem%elnod(e,n),i) =  fext_glob(elem%elnod(e,n),i) + elem%f_ext(e,n,i)
      end do !element row
    end do ! Element node
  end do ! e
end subroutine 



end module Matrices