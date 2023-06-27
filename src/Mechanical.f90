module Mechanical
use Domain
use ElementData
use NodeData
use Matrices

contains

!!!!!!!!!!!!!!!Gradv = L = dvx/dx dvx/dy  dvx/dz
!!!!!!!!!!!!!!!!!!!        dvy/dx dvy/dy  dvy/dz
!!!! E = 1/2 (L+LT)
!!!! R = 1/2 (L-LT)
!THIS SHOULD BE DONE AT t+1/2dt
subroutine cal_elem_strains ()
  implicit none
  integer :: e, i,j,k, gp, d, n
  real(fp_kind), dimension(dim,nodxelem) ::temp
  real(fp_kind) :: f
  
  elem%str_rate = 0.0d0
  elem%rot_rate = 0.0d0
  
  do e=1, elem_count
    do gp = 1, elem%gausspc(e)
      !Is only linear matrix?    
      !!!TODO: CHANGE FROM MATRIX OPERATION TO SIMPLE OPERATION
      f = 1.0d0/elem%detJ(e,gp)
      temp = elem%dHxy_detJ(e,gp,:,:) * f!!!!TODO: MODIFY BY MULTIPLYING
      elem%strain(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%uele (e,:,:)) 
      do n=1, nodxelem  
        do d=1, dim
          elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          elem%rot_rate(e,gp, d,d) = 0.0d0
        end do
        !!!! TO AVOID ALL MATMULT
        elem%str_rate(e,gp, 1,2) = elem%str_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) &!!!!dvx/dy
                                   + temp(1,n) * elem%vele (e,dim*(n-1)+2,1)
        elem%rot_rate(e,gp, 1,2) = elem%rot_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) & !!!!dvx/dx
                                   - temp(1,n) * elem%vele (e,dim*(n-1)+2,1)                           !!!!
        if (dim .eq. 3) then
          elem%str_rate(e,gp, 2,3) = elem%str_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &!!!d/dz*vy     
                                     + temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          elem%str_rate(e,gp, 1,3) = elem%str_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     + temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz     
          elem%rot_rate(e,gp, 2,3) = elem%rot_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &
                                     - temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          elem%rot_rate(e,gp, 1,3) = elem%str_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     - temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz    
        end if     
      end do !Nod x elem
      elem%str_rate(e,gp, 1,2) = 0.5 * elem%str_rate(e,gp, 1,2); 
      elem%rot_rate(e,gp, 1,2) = 0.5 * elem%rot_rate(e,gp, 1,2)      

      elem%str_rate(e,gp, 2,1) =     elem%str_rate(e,gp, 1,2)
      elem%rot_rate(e,gp, 2,1) =    -elem%rot_rate(e,gp, 1,2)
      if (dim .eq. 3) then
        elem%str_rate(e,gp, 1,3) = 0.5 * elem%str_rate(e,gp, 1,3); elem%str_rate(e,gp, 2,3) = 0.5 * elem%str_rate(e,gp, 2,3)
        elem%rot_rate(e,gp, 1,3) = 0.5 * elem%rot_rate(e,gp, 1,3); elem%rot_rate(e,gp, 2,3) = 0.5 * elem%rot_rate(e,gp, 2,3)
        
        elem%str_rate(e,gp, 3,2) =     elem%str_rate(e,gp, 2,3)
        elem%str_rate(e,gp, 3,1) =     elem%str_rate(e,gp, 1,3)

        elem%rot_rate(e,gp, 3,2) =     -elem%rot_rate(e,gp, 2,3)
        elem%rot_rate(e,gp, 3,1) =     -elem%rot_rate(e,gp, 1,3)
      end if

      !elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
      !print *, "strain rate ", elem%str_rate(e,gp,:,:)
      !print *, "rot    rate ", elem%rot_rate(e,gp,:,:)
    end do !gp
  end do !element
end subroutine

!calc int_forces
subroutine cal_elem_forces ()
  implicit none
  integer :: e, i,j,k, gp,n, d
  real(fp_kind), dimension(dim*nodxelem,1) ::f
  real(fp_kind) :: w
  elem%f_int = 0.0d0
  w = 1.0d0 !!! Full integration

  do e=1, elem_count
    if (elem%gausspc(e) .eq. 1) then
      w = 2.0d0**dim
    end if
    do gp = 1, elem%gausspc(e)
      !print *, "elem%dHxy_detJ(e,gp,1", elem%dHxy_detJ(e,gp,1,:)
      !print *, "elem%dHxy_detJ(e,gp,2", elem%dHxy_detJ(e,gp,2,:)
      do n=1, nodxelem
      !Is only linear matrix?    
      !elem%f_int(e,n,d) =  
      !f (:,:) = matmul(transpose(elem%bl(e,gp,:,:)),elem%sigma (e,:,:))
      !!!! TO AVOID MATRIX MULTIPLICATIONS (8x6 = 48 in bathe notation with several nonzeros)
      !!!!! F = BT x sigma = [dh1/dx dh1/dy ] x [ sxx sxy]
      !!!!!                = [dh2/dx dh2/dy ]   [ syx syy]
      !!!!! 
        do d=1, dim
          elem%f_int(e,n,d) = elem%f_int(e,n,d) + elem%dHxy_detJ(e,gp,d,n) * elem%sigma (e,gp, d,d)
        end do
        if (dim .eq. 2) then  !!!!! TODO: CHANGE WITH BENSON 1992 - EQ 2.4.2.11 FOR SIMPLICITY
          !!elem%f_int(e,n,1) = 
          elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) 
          elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2)
        else 
          elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) + &
                                                  elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 1,3)
          elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) + &
                                                  elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 2,3)
          elem%f_int(e,n,3) = elem%f_int(e,n,3) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 2,3) + &
                                                  elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,3)
        end if
        !print *, "Node ", n, "Force ", elem%f_int(e,n,:) 
      end do! nod x elem
    end do !gp
    elem%f_int(e,:,:) = elem%f_int(e,:,:) * w
  end do!elem
end subroutine


!!!! AFTER CALCULATING VELE 
!!!!!THIS HOURGLASS CONTROL IS BASED ON 
!!!!! GOUDREAU 1982 --> Used this simplified hourglass correction
!!!!! FLANAGAN 1981
!!!!! ATTENTION: IN FLANAGAN INTRINSIC COORDINATES ARE FROM -1/2 to 1/2
!!!!! SO: h1=1/4(1-2r)(1-2s) (Flanagan Eqn 55). 
!!!!! With our instrinsic from -1 to 1 , in  Table 2 Sigma should be 
!!!! Sigma is quadratic (2D) or qubic(3D) coefficient of 
subroutine calc_hourglass_forces
  implicit none
  integer :: e, n, j, d, gp, jmax
  real(fp_kind), dimension(4, nodxelem) :: Sig !! 4 COLUMNVECTORS IN 2D ONLY first is used
  real(fp_kind), dimension(nodxelem,dim):: vel!!!!DIFFERENT FROM vele which is an 8 x 1 vector
  real(fp_kind), dimension(dim,4) :: hmod
!real(fp_kind), dimension(1,4) :: test
  if (dim .eq. 2) then
    jmax = 1
  else
    jmax = 4
  end if
  !!! TODO: TO BE DEFINED IN DOMAIN ONCE 
  hmod(:,:) = 0.0d0
  elem%hourg_nodf(:,:,:) = 0.0d0
  if (dim .eq. 3) then
    Sig(1,:) = [ 0.125, 0.125,-0.125,-0.125,-0.125,-0.125, 0.125, 0.125]
    Sig(2,:) = [ 0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125,-0.125]
    Sig(3,:) = [ 0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125]
    Sig(4,:) = [-0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.125,-0.125]
  else 
    Sig(1,:) = [ 0.25, -0.25, 0.25,-0.25] !!!  
  end if
  
  gp = 1
  do e=1, elem_count    
    if (elem%gausspc(e) .eq. 1) then
      !test = matmul (elem%dHxy(e,gp,:,:),transpose(Sig(:,:))) !!!!SHOULD BE ORTHOGONAL
      !print *, "test ", test
      !print *, "dHxy ", elem%dHxy(e,gp,:,:)
          
      do n=1,nodxelem
        do d=1,dim
          vel (n,d) = nod%v(elem%elnod(e,n),d)    
        end do
      end do
      do j =1,jmax
        do n=1,nodxelem
          hmod(:,j) = hmod(:,j) + vel (n,:)*Sig(j,n) !!!!! ":" is on dimension, GOUDREAU EQN (30)
        end do
      end do
      
      do n=1,nodxelem
        do j =1,4  
            elem%hourg_nodf(e,n,:) = elem%hourg_nodf(e,n,:) - hmod(:,j)*Sig(j,n)
          end do
      end do
  else
    !print *, "NO HOURGLASS"
    end if
  end do !elemen
end subroutine

!!!!!!---------------------------
!!!!!!MIE-GRUNEISSEN EQN OF STATE
!!!!!!!--------------------------
!real(fp_kind) function EOS(EQ, Cs0, P00, Density, Density0
real(fp_kind) function EOS(Cs0, P00, Density, Density0)
  !integer, intent (in) :: EQ
  real(fp_kind), intent(in) :: Cs0, P00, Density, Density0
  ! switch (EQ)
    ! case 0:
      EOS =  P00+(Cs0*Cs0)*(Density-Density0)
      ! break;

    ! case 1:
      ! return P00+(Density0*Cs0*Cs0/7.0)*(pow(Density/Density0,7.0)-1);
      ! break;

    ! case 2:
      ! return (Cs0*Cs0)*Density;
      ! break;

    ! default:
      ! std::cout << "Please correct Pressure Equation No and run again" << std::endl;
      ! std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
      ! std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
      ! std::cout << "2 => (Cs*Cs)*Density" << std::endl;
      ! abort();
      ! break;
  
end function EOS

!!!!! ASSUME VOLUME IS ALREADY CALCULATED
subroutine calc_elem_density ()
  implicit none
  real(fp_kind), dimension(dim,dim) :: F !def gradient
  real(fp_kind), dimension(nodxelem,dim) :: x !def gradient
  
  integer :: e, n, gp
  do e = 1, elem_count

    if (elem%gausspc(e) .eq. 1) then
      elem%rho(e,1) = elem%rho_0(e,1)*elem%vol_0(e)/elem%vol(e) !IS THE SAME
    else 
      do n=1,nodxelem 
        x(n,:) = nod%x(elem%elnod(e,n),:) !!!CURRENT CONFIG
      end do
      do gp = 1, elem%gausspc(e)
        !!!!CALCULATE DEFORMATION GRADIENT
        F = matmul(elem%dHxy0(e,gp,:,:),x) !!!!TODO; DO THIS FOR ALL 
        !print *, "det F", det(F)
        !elem%rho(e,gp) = elem%rho_0(e,gp)* elem%detJ(e,gp)
        elem%rho(e,gp) = elem%rho_0(e,gp)* det(F)
      end do
    end if
  end do
end subroutine

subroutine calc_elem_pressure ()
  implicit none
  integer :: e, gp
  
  gp = 1
  do e = 1, elem_count
    do gp = 1, elem%gausspc(e)
      !print *, "cs rho rho0 ", elem%cs(e), elem%rho(e,gp),elem%rho_0(e,gp)
      elem%pressure(e,gp) = EOS(elem%cs(e),0.0d0,elem%rho(e,gp),elem%rho_0(e,gp))
    end do
  end do
end subroutine


!!!!!! IT ASSUMES PRESSURE AND STRAIN RATES ARE ALREADY CALCULATED
!!!!!! (AT t+1/2 to avoid stress at rigid rotations, see Benson 1992)
subroutine CalcStressStrain (dt) 

  implicit none
  real(fp_kind) :: SRT(3,3), RS(3,3), ident(3,3)
  integer :: e,gp
  real(fp_kind) ,intent(in):: dt
  
  real(fp_kind) :: p00
  
  p00 = 0.
  
  ident = 0.0d0
  ident (1,1) = 1.0d0; ident (2,2) = 1.0d0; ident (3,3) = 1.0d0
  
  ! !!!$omp parallel do num_threads(Nproc) private (RotationRateT, Stress, SRT, RS)
  do e = 1, elem_count
    
    do gp=1,elem%gausspc(e)
    !!!!! ALLOCATE REDUCED VECTOR TO TENSOR 
    ! do d=1,dim
! !      stress(i,i) = elem%sigma (e,gp,i
      ! str_rate(d,d) = elem%str_rate(e,gp,d,1)
      ! rot_rate(d,d) = elem%rot_rate(e,gp,d,1)
    ! end do
    ! str_rate(1,2) = elem%str_rate(e,gp,dim+1,1)/2.0 ; str_rate(2,1) = str_rate(1,2)
    ! rot_rate(1,2) = elem%rot_rate(e,gp,dim+1,1)/2.0 ; rot_rate(2,1) = rot_rate(1,2)
    ! if (dim .eq. 3) then
      ! str_rate(2,3) = elem%str_rate(e,gp,dim+2,1)/2.0 ; str_rate(3,2) = str_rate(2,3)   !!!BATHE table 6.6 p. 556
      ! str_rate(3,1) = elem%str_rate(e,gp,dim+3,1)/2.0 ; str_rate(1,3) = str_rate(3,1)

      ! ! str_rate(2,3) = elem%str_rate(e,gp,dim+2,1) ; str_rate(3,2) = str_rate(2,3)   !!!BATHE table 6.6 p. 556
      ! ! str_rate(3,1) = elem%str_rate(e,gp,dim+3,1) ; str_rate(1,3) = str_rate(3,1)
    ! end if
    ! pt%pressure(i) = EOS(0, pt%cs(i), p00,pt%rho(i), pt%rho_0(i))
    ! if (i==52) then
    ! !print *, "pt%pressure(i)", pt%pressure(i),", cs ", pt%cs(i), "p00", p00, ", rho", p00,pt%rho(i), ", rho 0", p00,pt%rho_0(i)
    ! end if
    ! RotationRateT = transpose (elem%rot_rate(e,:,:))

      SRT = MatMul(elem%shear_stress(e,gp,:,:),transpose(elem%rot_rate(e,gp,:,:)))
      RS  = MatMul(elem%rot_rate(e,gp,:,:), elem%shear_stress(e,gp,:,:))
      
      ! !print *, "RS", RS
      !print *, "mat g", mat_G
      elem%shear_stress(e,gp,:,:)	= dt * (2.0 * mat_G *(elem%str_rate(e,gp,:,:) - 1.0/3.0 * &
                                   (elem%str_rate(e,gp,1,1)+elem%str_rate(e,gp,2,2)+elem%str_rate(e,gp,3,3))*ident) &
                                   +SRT+RS) + elem%shear_stress(e,gp,:,:)
      elem%sigma(e,gp,:,:) = -elem%pressure(e,gp) * ident + elem%shear_stress(e,gp,:,:)	!Fraser, eq 3.32
      !print *, "elem ", e, ", sigma ", elem%sigma(e,gp,:,:)
    ! !pt%strain(i)			= dt*pt%str_rate(i + Strain;
    end do !gauss point
  end do
  ! !!!!$omp end parallel do    


end subroutine CalcStressStrain

subroutine calc_elem_vol ()
  implicit none
  integer :: e, gp
  real(fp_kind):: w
  
  ! P00+(Cs0*Cs0)*(Density-Density0);
  do e = 1, elem_count
    elem%vol(e) = 0.0d0
    !print *, "elem%gausspc(e)", elem%gausspc(e)
    if (elem%gausspc(e).eq.1) then
      w = 2.0**dim
    else if (elem%gausspc(e).eq. 8 ) then
      w = 1.0
    end if
    do gp=1,elem%gausspc(e)
      !elem%vol(e) = 
      print *, "elem e j, w", elem%detJ(e,gp), w
      elem%vol(e) = elem%vol(e) + elem%detJ(e,gp)*w
    end do !gp  
    !print *, "Elem ", e, "vol ",elem%vol(e)
  end do

end subroutine




subroutine impose_bcv
  implicit none
  integer :: n, d
  n = 1
  do while (n .le. node_count)    
    do d=1,dim
      if (nod%is_bcv(n,d) .eqv. .true.) then
        nod%v(n,d) = nod%bcv(n,d)
        !print *, "nod ", n, ", ",nod%bcv(n,d), ", d", d
      end if
      
      if (nod%is_fix(n,d) .eqv. .true.) then
        nod%v(n,d) = 0.0d0
      end if 

    end do !dim
    n = n + 1
  end do !Node    
end subroutine

subroutine impose_bca
  implicit none
  integer :: n, d
  do n=1,node_count    
    do d =1, dim
      ! if (nod%is_bcv(n,d) .eqv. .true.) then
        ! nod%v(n,d) = nod%bcv(n,d)
        ! print *, "nod ", n, ", ",nod%bcv(n,d), ", d", d
      ! end if      
      if (nod%is_fix(n,d) .eqv. .true.) then
        nod%a(n,d) = 0.0d0
      end if 
    end do !dim
  end do !Node    
end subroutine

end module
