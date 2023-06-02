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
    elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
  end do
end subroutine

!calc int_forces
subroutine cal_elem_forces ()
  implicit none
  integer :: e, i,j,k, gp,n, d
  real(fp_kind), dimension(dim*nodxelem,1) ::f
  
  gp = 1
  do e=1, elem_count
    do n=1, nodxelem
    !Is only linear matrix?    
    !elem%f_int(e,n,d) =  
    !f (:,:) = matmul(transpose(elem%bl(e,gp,:,:)),elem%sigma (e,:,:))
    !!!! TO AVOID MATRIX MULTIPLICATIONS (8x6 = 48 in bathe notation with several nonzeros)
      do d=1, dim
        elem%f_int(e,n,d) = elem%dHxy(e,gp,n,d) * elem%tau (e,gp, d,d)
      end do
      elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy(e,gp,2,n) * elem%tau (e,gp, 1,2) + &
                                              elem%dHxy(e,gp,3,n) * elem%tau (e,gp, 1,3)
      elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy(e,gp,1,n) * elem%tau (e,gp, 1,2) + &
                                              elem%dHxy(e,gp,3,n) * elem%tau (e,gp, 2,3)
      elem%f_int(e,n,3) = elem%f_int(e,n,3) + elem%dHxy(e,gp,2,n) * elem%tau (e,gp, 2,3) + &
                                              elem%dHxy(e,gp,1,n) * elem%tau (e,gp, 1,3)
    end do
  end do
end subroutine


!!!! AFTER CALCULATING VELE 
subroutine calc_hourglass_forces
  implicit none
  integer :: e, n, j, d, gp
  real(fp_kind), dimension(4, nodxelem) :: Sig !! 4 COLUMNVECTORS IN 2D ONLY first is used
  real(fp_kind), dimension(nodxelem,dim):: vel!!!!DIFFERENT FROM vele which is an 8 x 1 vector
  real(fp_kind), dimension(dim,4) :: hmod
!real(fp_kind), dimension(1,4) :: test
  
  
  if (dim .eq. 3) then
    Sig(1,:) = [ 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0]
    Sig(2,:) = [ 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0]
    Sig(3,:) = [ 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0]
    Sig(4,:) = [-1.0, 1.0,-1.0, 1.0, 1.0,-1.0, 1.0,-1.0]
  end if
  
  gp = 1
  do e=1, elem_count    
    !test = matmul (elem%dHxy(e,gp,:,:),transpose(Sig(:,:))) !!!!SHOULD BE ORTHOGONAL
    !print *, "test ", test
    !print *, "dHxy ", elem%dHxy(e,gp,:,:)
        
    do n=1,nodxelem
      do d=1,dim
        vel (n,d) = nod%v(elem%elnod(e,n),d)    
      end do
    end do
    do j =1,4
      do n=1,nodxelem
        hmod(:,j) = hmod(:,j) + vel (n,:)*Sig(j,n) !!!!! ":" is on dimension
      end do
    end do
    
    do n=1,nodxelem
      do j =1,4  
          elem%hourg_nodf(e,n,:) = elem%hourg_nodf(e,n,:) - hmod(:,j)*Sig(j,n)
        end do
    end do
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
  integer :: e
  do e = 1, elem_count
    elem%rho(e) = elem%mass(e)/elem%vol(e)
  end do
end subroutine

subroutine calc_elem_pressure ()
  implicit none
  integer :: e, gp
  
  gp = 1
  do e = 1, elem_count
    elem%pressure(e,gp) = EOS(elem%cs(e),0.0d0,elem%rho(e),elem%rho_0(e,gp))
  end do
end subroutine

subroutine CalcStressStrain (dt) 

  implicit none
  real(fp_kind) :: RotationRateT(3,3), Stress(3,3), SRT(3,3), RS(3,3), ident(3,3)
  integer :: i
  real(fp_kind) ,intent(in):: dt
  
  real(fp_kind) :: p00
  
  p00 = 0.
  
  ident = 0.0d0
  ident (1,1) = 1.0d0; ident (2,2) = 1.0d0; ident (3,3) = 1.0d0
  
  ! !!!$omp parallel do num_threads(Nproc) private (RotationRateT, Stress, SRT, RS)
  do i = 1, elem_count
    ! pt%pressure(i) = EOS(0, pt%cs(i), p00,pt%rho(i), pt%rho_0(i))
    ! if (i==52) then
    ! !print *, "pt%pressure(i)", pt%pressure(i),", cs ", pt%cs(i), "p00", p00, ", rho", p00,pt%rho(i), ", rho 0", p00,pt%rho_0(i)
    ! end if
    ! RotationRateT = transpose (elem%rot_rate(e,:,:))

    ! SRT = MatMul(pt%shear_stress(i,:,:),RotationRateT)
    ! RS  = MatMul(pt%rot_rate(i,:,:), pt%shear_stress(i,:,:))
    
    ! !print *, "RS", RS
    ! pt%shear_stress(i,:,:)	= dt * (2.0 * mat_G *(pt%str_rate(i,:,:)-1.0/3.0 * &
                                 ! (pt%str_rate(i,1,1)+pt%str_rate(i,2,2)+pt%str_rate(i,3,3))*ident) &
                                 ! +SRT+RS) + pt%shear_stress(i,:,:)
    ! pt%sigma(i,:,:)			= -pt%pressure(i) * ident + pt%shear_stress(i,:,:)	!Fraser, eq 3.32
    ! !print *, "particle ", i, ", rot_rate ", pt%rot_rate(i,:,:)
    ! !pt%strain(i)			= dt*pt%str_rate(i + Strain;
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
    if (elem%gausspc(e).eq.1) then
      w = 8.0
    end if
    do gp=1,elem%gausspc(e)
      !elem%vol(e) = 
      print *, "elem e j", elem%detJ(e,gp)
      elem%vol(e) = elem%vol(e) + elem%detJ(e,gp)*w
      end do !gp  
    print *, "Elem ", e, "vol ",elem%vol(e)
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