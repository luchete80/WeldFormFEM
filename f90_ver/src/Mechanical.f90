module Mechanical
use Domain
use ElementData
use NodeData
use Matrices
use mymath  !TRACE

contains

!!!!!!!!!!!!!!!Gradv = L = dvx/dx dvx/dy  dvx/dz
!!!!!!!!!!!!!!!!!!!        dvy/dx dvy/dy  dvy/dz
!!!! E = 1/2 (L+LT)
!!!! R = 1/2 (L-LT)
!THIS SHOULD BE DONE AT t+1/2dt
subroutine cal_elem_strains ()
  use omp_lib
  
  implicit none
  integer :: e, i,j,k, gp, d, n
  real(fp_kind), dimension(dim,nodxelem) ::temp
  real(fp_kind) :: f
  real(fp_kind) :: test(1,6),test33(3,3) !ifwanted to test in tensor form
  
  elem%str_rate = 0.0d0
  elem%rot_rate = 0.0d0
  
  !$omp parallel do num_threads(Nproc) private (e,gp,d, temp,n) 
  do e=1, elem_count
    do gp = 1, elem%gausspc(e)
      !Is only linear matrix?    
      !!!TODO: CHANGE FROM MATRIX OPERATION TO SIMPLE OPERATION
      f = 1.0d0/elem%detJ(e,gp)
      temp = elem%dHxy_detJ(e,gp,:,:) * f!!!!TODO: MODIFY BY MULTIPLYING
      ! elem%strain(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%uele (e,:,:)) 
      ! !print *, "standard stran rate calc (matricial) "
      ! !!!!!! DEFAULT (TODO: CHECK IF IS SLOW)
      ! test = f* matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:))  ! (6x24)(24x1)
      ! !!print *, "e11 e22 e33 2e12 2e23 2e31", test
      ! test33(1,1) = test(1,1);test33(2,2) = test(1,2);test33(3,3) = test(1,3);
      ! test33(1,2) = test(1,4)*0.5;test33(2,1) =test33(1,2);
      ! test33(2,3) = test(1,5)*0.5;test33(3,2) =test33(2,3);
      ! test33(3,1) = test(1,6)*0.5;test33(1,3) =test33(3,1);
      

      ! test33 = 0.5*(test33+transpose(test33));
      ! ! !print *, "str rate test", test33
      
      ! ! test33 = 0.5*(test33-transpose(test33));
      ! ! print *, "rot rate test", test33
      !!!!! er = B x U = [dh1/dx  0      dh1/dy  0 ] x [ U1 v1 u2 ]T
      !!!!!              [0       dh1/dy 0      dh2/dy ]   
      !!!!!              [dh1/dy  dh1/dx dh2/dy  dh2/dx ]                     
      do n=1, nodxelem  
        do d=1, dim
          !print *, "node dim dHxy vele", n,d,temp(d,n) , elem%vele (e,dim*(n-1)+d,1) 
          elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          elem%rot_rate(e,gp, d,d) = 0.0d0
        end do
        !!!! TO AVOID ALL MATMULT
        elem%str_rate(e,gp, 1,2) = elem%str_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) &!!!!dvx/dy
                                   + temp(1,n) * elem%vele (e,dim*(n-1)+2,1)
        elem%rot_rate(e,gp, 1,2) = elem%rot_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) & !!!!dvx/dy
                                   - temp(1,n) * elem%vele (e,dim*(n-1)+2,1)                           !!!!
        
        !!! er hoop = vr/r
        if (dim .eq. 2 .and. bind_dom_type .eq. 3) then 
          elem%str_rate(e,gp, 3,3) = elem%str_rate(e,gp, 3,3) + elem%vele (e,dim*(n-1)+1,1) / elem%radius(e,gp)
          elem%rot_rate(e,gp, 3,3) = 0.0d0
        end if 
        
        if (dim == 3) then
          elem%str_rate(e,gp, 2,3) = elem%str_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &!!!d/dz*vy     
                                     + temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          elem%str_rate(e,gp, 1,3) = elem%str_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     + temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz     
          elem%rot_rate(e,gp, 2,3) = elem%rot_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &
                                     - temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          elem%rot_rate(e,gp, 1,3) = elem%rot_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
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
			
			!!! IF COMPLETE MULTIPLICATION (MORE CALC)
      !!!elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
      !print *, "simlpified strain rate "
      !print *, "strain rate ", elem%str_rate(e,gp,:,:)
      !print *, "rot    rate ", elem%rot_rate(e,gp,:,:)
    end do !gp
  end do !element
  !$omp end parallel do    
end subroutine

!!!!! TWO WAYS:FROM STRAIN RATE OR BEING U FROM: F = U R 
subroutine cal_elem_strain_inc (dt)
  real(fp_kind), intent(in) :: dt
  do e=1, elem_count
    do gp = 1, elem%gausspc(e)
      elem%str_inc (e,gp, :,:) = elem%str_rate(e,gp,:,:) * dt
    end do !gp
  end do !element
end subroutine

!calc int_forces
subroutine cal_elem_forces ()
	use omp_lib
	
  implicit none
  integer :: e, i,j,k, gp,n, d
  real(fp_kind), dimension(dim*nodxelem,1) ::f
  real(fp_kind) :: w
  !TESTING
  real (fp_kind) :: sigma_test(6,1) !ORDERED
  real(fp_kind) :: test(24,1) !ifwanted to test in tensor form
  real(fp_kind) :: area, f2 ! Axisymm
  
  f2 = 1.0d0 !!!! AXISYMM FACTOR IN CASE OF INTGRAL CASE
  
  elem%f_int = 0.0d0
  w = 1.0d0 !!! Full integration
	
	! !$omp parallel do num_threads(Nproc) private (e,gp,d, w,n) 
  do e=1, elem_count
    if (elem%gausspc(e) .eq. 1) then
      w = 2.0d0**dim
    end if
    do gp = 1, elem%gausspc(e)
      !print *, "elem%dHxy_detJ(e,gp,1", elem%dHxy_detJ(e,gp,1,:)
      !print *, "elem%dHxy_detJ(e,gp,2", elem%dHxy_detJ(e,gp,2,:)
      ! sigma_test (:,1)=[elem%sigma (e,gp, 1,1),elem%sigma (e,gp, 2,2),elem%sigma (e,gp, 3,3),&
                        ! elem%sigma (e,gp, 1,2),elem%sigma (e,gp, 2,3),elem%sigma (e,gp, 3,1)]
      ! test = w*matmul(transpose(elem%bl(e,gp,:,:)),sigma_test)  ! (24x6)(6x1)
      !print *, "test force", test
      
      !print *, "dHdxy, 1", elem%dHxy_detJ(e,gp,1,:)
      !print *, "dHdxy, 2", elem%dHxy_detJ(e,gp,2,:)
      !print *, "dHdxy, 3", elem%dHxy_detJ(e,gp,1,:)
      
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
          if (bind_dom_type .eq. 3 .and. axisymm_vol_weight .eqv. .true.) then
            !print *, "AAAAAAAAAAAAAAAAAAA", elem%radius(e,gp)
            f2 = elem%radius(e,gp)
          end if
          !!elem%f_int(e,n,1) = 
          elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) * f2
          elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) * f2
          
          !!!These are dividing per r so in the vol weight is like this
          !!! Goudreau 1982 eq. 19
           
          if (bind_dom_type .eq. 3 .and. axisymm_vol_weight .eqv. .true.) then
            elem%f_int(e,n,1) = elem%f_int(e,n,1) + (elem%sigma (e,gp, 1,1) - &
                                                     elem%sigma (e,gp, 3,3) ) * elem%detJ(e,gp)
            elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%sigma (e,gp, 1,2) * elem%detJ(e,gp)
          end if
        else 
          elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) + &
                                                  elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 1,3)
          elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) + &
                                                  elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 2,3)
          elem%f_int(e,n,3) = elem%f_int(e,n,3) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 2,3) + &
                                                  elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,3)
        end if
        !print *, "Element force Node ", n, "F  ", elem%f_int(e,n,:) * w 
      end do! nod x elem
      !print *, "test ", w * elem%dHxy_detJ(e,gp,3,8)  * elem%sigma (e,gp, 3,3)
      !print *, "dHxy ", elem%dHxy_detJ(e,gp,3,8), "w ", w
      !print *, "s33 ", elem%sigma (e,gp, 3,3)
    end do !gp
    elem%f_int(e,:,:) = elem%f_int(e,:,:) * w
  end do!elem
	! !$omp end parallel do    
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
  real(fp_kind) :: c_h
  real(fp_kind), dimension(4, nodxelem) :: Sig !! 4 COLUMNVECTORS IN 2D ONLY first is used
  real(fp_kind), dimension(nodxelem,dim):: vel!!!!DIFFERENT FROM vele which is an 8 x 1 vector
  real(fp_kind), dimension(dim,4) :: hmod
  real(fp_kind) :: test
!real(fp_kind), dimension(1,4) :: test
  if (dim .eq. 2) then
    jmax = 1
  else
    jmax = 4
  end if
  !!! TODO: TO BE DEFINED IN DOMAIN ONCE hmod(:,:) = 0.0d0
  elem%hourg_nodf(:,:,:) = 0.0d0
  if (dim .eq. 3) then
    !Also in Flanagan Appendix (data GB in program)
    Sig(1,:) = [ 0.125, 0.125,-0.125,-0.125,-0.125,-0.125, 0.125, 0.125]
    Sig(2,:) = [ 0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125,-0.125]
    Sig(3,:) = [ 0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125]
    Sig(4,:) = [-0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.125,-0.125]
    Sig(:,:) = Sig(:,:) * 8
  else 
    Sig(1,:) = [ 0.25, -0.25, 0.25,-0.25] !!!  
    Sig(:,:) = Sig(:,:) * 4
  end if
  
  gp = 1
  do e=1, elem_count    
    if (elem%gausspc(e) .eq. 1) then
      hmod(:,:) = 0.0d0
      !test = matmul (elem%dHxy(e,gp,:,:),transpose(Sig(:,:))) !!!!SHOULD BE ORTHOGONAL
      !print *, "test ", test
      !print *, "dHxy ", elem%dHxy(e,gp,:,:)
          
      do n=1,nodxelem
        do d=1,dim
          vel (n,d) = nod%v(elem%elnod(e,n),d)    
        end do
      end do
      do j =1,jmax !MODE
        do n=1,nodxelem !1:4 or 8
          !print *, "elem v ", vel (n,:)
          hmod(:,j) = hmod(:,j) + vel (n,:)*Sig(j,n) !!!!! ":" is on dimension d, GOUDREAU EQN (30)
        end do
      end do
      
      
      
      !!!!!!!!! GOUDREAU 1982
      do n=1,nodxelem
        do j =1,jmax 
            elem%hourg_nodf(e,n,:) = elem%hourg_nodf(e,n,:) - hmod(:,j)*Sig(j,n)
          end do
      end do
      c_h  = 0.06 * elem%vol(e)**(0.6666666) * elem%rho(e,1) * 0.25 * mat_cs0
      
      !print *, "hourglass c ", c_h
      elem%hourg_nodf(e,:,:) = elem%hourg_nodf(e,:,:) * c_h
      
      
      !!!!!!!!! FLANAGAN 1981
      
      
      ! do n=1,nodxelem
      ! print *, "hourglass forces", elem%hourg_nodf(e,n,:) 
      ! end do
      ! ! elem%hourg_nodf(e,:,:) = - matmul(matmul(transpose(Sig(:,:)),Sig(:,:)),vel (:,:)) * c_h 
      
      ! print *, "alt hourglass forces", elem%hourg_nodf(e,:,:) 
      ! do n=1,nodxelem
        ! test = 0
        ! do d=1,3
          ! test = test + elem%hourg_nodf(e,n,d) *  vel (n,d)
        ! end do
        ! print *, "dot vel and hg forces ", test
      ! end do

      ! !dimension and modes
      ! do d=1,3
        ! do j =1,jmax !MODE
          ! test = 0
          ! do n=1,nodxelem !1:4 or 8
            ! !print *, "elem v ", vel (n,:)
            ! !hmod(:,j) = hmod(:,j) + vel (n,:)*Sig(j,n) !!!!! ":" is on dimension d, GOUDREAU EQN (30)
            ! test = test + vel (n,d)*Sig(j,n) !!!!! ":" is on dimension d, GOUDREAU EQN (30)
          ! end do
          ! ! print *, "Mode ", j, "Dim", d
          ! print *, "mode", j, "dim ", d, "sum v x Sig", test, ", hmod ", hmod(d,j)
          ! ! print *, "force "
        ! end do
      ! end do

      ! do n=1,nodxelem
          ! print *, "v ", vel (n,:)
      ! end do
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
    do gp=1, elem%gausspc(e)
    !if (elem%gausspc(e) .eq. 1) then
      elem%rho(e,gp) = elem%rho_0(e,gp)*elem%vol_0(e)/elem%vol(e) !IS THE SAME
    end do
  end do
end subroutine


!!!!! ONLY FOR GREEN NAGHDI OR COMPUTE JAUMANN FROM DL = I -FINV
!!!!! JAUMANN IS BETTER TO BE CALCULATED USING W FROM STRAIN RATE
!!!!! FIRST SHOULD BE CALCULATED DEFORMATION GRADIENT 
subroutine calc_def_grad ()
  implicit none
  real(fp_kind), dimension(dim,dim) :: F,ident(3,3) !def gradient
  real(fp_kind), dimension(nodxelem,dim) :: x,u !def gradient
  integer :: e, n, gp
  
  ident = 0.0d0
  ident (1,1) = 1.0d0; ident (2,2) = 1.0d0; ident (3,3) = 1.0d0
  
  do e = 1, elem_count
    do n=1,nodxelem 
      !x(n,:) = nod%x(elem%elnod(e,n),:) !!!CURRENT CONFIG
      !u(n,:) = nod%u(elem%elnod(e,n),:) !!!CURRENT CONFIG
      u(n,:) = nod%u_inc(elem%elnod(e,n),:) !!!CURRENT CONFIG
    end do
    do gp = 1, elem%gausspc(e)
      !F = matmul(elem%dHxy0(e,gp,:,:),x) !!!!TODO; DO THIS FOR ALL 
      F = ident + 1.0d0/elem%detJ(e,gp) * matmul(elem%dHxy_detJ(e,gp,:,:),u) !!!!TODO; DO THIS FOR ALL 
      elem%def_grad(e,gp,:,:) = F
      if (gp==1) then
      ! print *, "def grad F", F(1,1), F(1,2), F(1,3)
      ! print *, F(2,1), F(2,2), F(2,3)
      ! print *, F(3,1), F(3,2), F(3,3)
      end if
    end do

  end do!elem
end subroutine

subroutine calc_polar_urmat 
  real(fp_kind) ::U(6)
  do e = 1, elem_count
    do gp = 1, elem%gausspc(e)
      call polarCuppen(elem%def_grad(e,gp,:,:),U,elem%rmat(e,gp,:,:))
      elem%umat(e,gp,1,1) = U(1);elem%umat(e,gp,1,2) = U(2);elem%umat(e,gp,1,3) = U(3);
      elem%umat(e,gp,2,2) = U(4);elem%umat(e,gp,2,3) = U(5);elem%umat(e,gp,3,3) = U(6);
      elem%umat(e,gp,2,1) = elem%umat(e,gp,1,2);elem%umat(e,gp,3,1) = elem%umat(e,gp,1,3);elem%umat(e,gp,3,2) = elem%umat(e,gp,2,3); 
      if (gp == 1) then
        ! print *, "U ", U
        ! print *, "R ", elem%rmat(e,gp,:,:)
      end if 
    end do 
  end do
endsubroutine

subroutine cal_elem_strain_inc_from_umat
  do e=1, elem_count
    do gp = 1, elem%gausspc(e)
      elem%str_inc (e,gp, :,:)= elem%umat(e,gp,:,:)
      ! print *, "elem%str_inc (e,gp, :,:)", elem%str_inc (e,gp, :,:)
    end do !gp
  end do !element
end subroutine

subroutine cal_elem_strain_inc_from_str_rate (dt)
  use omp_lib
	real(fp_kind) :: dt
	integer::e
  !$omp parallel do num_threads(Nproc) private (e,gp)
	do e=1, elem_count
    do gp = 1, elem%gausspc(e)
      elem%str_inc (e,gp, :,:)= elem%str_rate(e,gp,:,:) * dt
      ! print *, "elem%str_inc (e,gp, :,:)", elem%str_inc (e,gp, :,:)
    end do !gp
  end do !element
	!$omp end parallel do
end subroutine

subroutine calc_inv_def_grad
  real(fp_kind) :: inv_f_incr(3,3)
  do e = 1, elem_count
    do gp = 1, elem%gausspc(e)
      !inv_f_incr = ident + 1.0d0/elem%detJ(e,gp) * matmul(elem%dHxy_detJ(e,gp,:,:),nod%x_prev()) !!!!TODO; DO THIS FOR ALL    
    end do
  end do
end subroutine calc_inv_def_grad

subroutine calc_strain_inc
  do e = 1, elem_count
    do gp = 1, elem%gausspc(e)
      
    end do
  end do
end subroutine calc_strain_inc

!!!!! ONLY FOR GREEN-NAHGDI
subroutine polardecomp()
  ! real(fp_kind), dimension(dim,dim) :: F
  ! real(fp_kind), Fsym(6)
  ! do e=1, elem_count
    ! do gp = 1, elem%gausspc(e)
      ! F = elem%def_grad(e,gp,:,:)
      ! Fsym(1) = F(1,1); Fsym(2) = F(1,2); Fsym(3) = F(1,3); 
      ! Fsym(4) = F(2,2); Fsym(5) = F(2,3); Fsym(6) = F(3,3); 
      ! !call polarCuppen(F, U, R)
    ! end do 
  ! end do
end subroutine


subroutine calc_elem_pressure ()
  call calc_elem_pressure_EOS
end subroutine

subroutine calc_elem_pressure_EOS ()
  implicit none
  integer :: e, gp
  
  gp = 1
  do e = 1, elem_count
    do gp = 1, elem%gausspc(e)
      ! print *, "cs rho rho0 ", elem%cs(e), elem%rho(e,gp),elem%rho_0(e,gp)
      elem%pressure(e,gp) = EOS(elem%cs(e),0.0d0,elem%rho(e,gp),elem%rho_0(e,gp))
    end do
  end do
end subroutine

!!!! STANDARD FORM- NON HYDRODYNAMIC WITH EOS
  ! double pressureIncrement = 0.0;

  ! for (intPointId = 0;
    ! pressureIncrement += getIntegrationPoint(intPointId)->StrainInc.trace();
  ! pressureIncrement /= getNumberOfIntegrationPoints();

  ! for (intPointId = 0
    ! getIntegrationPoint(intPointId)->pressure = getIntegrationPoint(intPointId)->Stress.thirdTrace() + K * pressureIncrement;
!    // Polar decomposition
!    F.polarCuppenLnU(_integrationPoint->StrainInc, _integrationPoint->R);    
! Assumes const module K (TODO: pass to element pointer)
! Needs to be calculated AFTER CALC STRAIN_INC FROM STRAIN RATE
subroutine calc_elem_pressure_from_strain (modK)
  implicit none
  real(fp_kind), intent(in) :: modK
  integer :: e, gp  
  real(fp_kind) :: press_inc
  
  gp = 1
  do e = 1, elem_count
    press_inc = 0.0d0
    do gp = 1, elem%gausspc(e)
      press_inc = press_inc + trace (elem%str_inc(e,gp,:,:))
    end do
    press_inc = -press_inc/elem%gausspc(e)
    ! print *, "press inc ", press_inc
    do gp = 1, elem%gausspc(e)    
          elem%pressure(e,gp) = -1.0/3.0 * trace (elem%sigma(e,gp,:,:)) + press_inc * modK
    end do
    ! print *, "mod K", modK
    ! print *, "strain inc", elem%str_inc(e,gp,:,:)
    ! print *, "press_inc ", press_inc
    ! print *, "elem%pressure(e,gp) FROM STRAIN", elem%pressure(e,1)
  end do
end subroutine

!TODO: MMAKE A MATERIAL OBJECT
!ASUMES constant K
subroutine calc_elem_wave_speed (modK)
  implicit none
  integer :: e, gp
  real(fp_kind), intent(in) :: modK
  do e = 1, elem_count
    do gp = 1, elem%gausspc(e)
      elem%c_s(e,gp) = sqrt(modK/elem%rho(e,gp))
       print *, " elem cs  ", elem%c_s(e,gp)
    end do
  end do  
end subroutine

!!!!!ACCORDING TO ABAQUS
!!!!!! p v1 = b1 rho cd Le e_vol_dot
!!!!!! p v2 = rho (b2 Le e_vol_dot)**2
!!!!!! SECOND ONE IS ONLY CALCULATED COMPRESSIVE LOADS
!!!!! AND ACCORDING TO BENSON 1992

subroutine calc_elem_shock_visc (dt)
  implicit none
  integer :: e, gp
  real(fp_kind), intent(in) :: dt
  real(fp_kind) :: vdot
  
  real(fp_kind) :: press_inc
  gp = 1
  do e = 1, elem_count
    do gp = 1, elem%gausspc(e)
      vdot = elem%vol_inc(e) / dt
      elem%p_visc(e,gp) = 0.06 * elem%rho(e,gp) * elem%c_s(e,gp) * vdot + 1.5 * elem%rho(e,gp) * vdot ** 2.0
      print *, "elem%p_visc(e,gp) tot", elem%p_visc(e,gp), "linear ", 1.5 * elem%rho(e,gp) * vdot ** 2.0
      print *, " elem cs  ", elem%c_s(e,gp)
    end do
  end do
end subroutine

subroutine calc_time_step
  ! // initialisation de la longueur caracteristique de l'element 0
  ! characteristicLength = elements(0)->getCharacteristicLength();

  ! // vitesse du son dans l'element 0
  ! elongationWaveSpeed = elements(0)->getElongationWaveSpeed();

  ! // valeur critique du timeStep step
  ! criticalTimeStep = characteristicLength / elongationWaveSpeed;

  ! for (long elementId = 1; elementId < elements.size(); elementId++)
  ! {
    ! // longueur caracteristique de l'element
    ! characteristicLength = elements(elementId)->getCharacteristicLength();

    ! // vitesse du son dans l'element
    ! elongationWaveSpeed = elements(elementId)->getElongationWaveSpeed();

    ! // valeur critique
    ! timeStep = characteristicLength / elongationWaveSpeed;

    ! // minimum des valeurs
    ! if (timeStep < criticalTimeStep)
      ! criticalTimeStep = timeStep;
  ! }
  ! return criticalTimeStep;

end subroutine


! //-----------------------------------------------------------------------------
! double ElHex8N3D::getCharacteristicLength()
! //-----------------------------------------------------------------------------
! {
    ! double A1, A2, A3, A4, A5, A6;
    ! AERT(A1, 0, 1, 2, 3)
    ! AERT(A2, 0, 1, 5, 4)
    ! AERT(A3, 1, 2, 6, 5)
    ! AERT(A4, 3, 2, 6, 7)
    ! AERT(A5, 0, 3, 7, 4)
    ! AERT(A6, 4, 5, 6, 7)
    ! A1 = dnlMax3(A1, A2, A3);
    ! A2 = dnlMax3(A4, A5, A6);
    ! A1 = dnlMax(A1, A2);
    ! return getVolume() / A1;
! }

!!!! ACORDING TO STANDARD SOLID FEM, PRESSURE FROM strain rate increment
! //-----------------------------------------------------------------------------
! void Element::computePressure()
! //-----------------------------------------------------------------------------
! {
  ! short intPointId;
  ! double K = material->getBulkModulus();
  ! double pressureIncrement = 0.0;

  ! for (intPointId = 0; intPointId < getNumberOfIntegrationPoints(); intPointId++)
  ! {
    ! pressureIncrement += getIntegrationPoint(intPointId)->StrainInc.trace();
  ! }

  ! pressureIncrement /= getNumberOfIntegrationPoints();

  ! for (intPointId = 0; intPointId < getNumberOfIntegrationPoints(); intPointId++)
  ! {
    ! getIntegrationPoint(intPointId)->pressure = getIntegrationPoint(intPointId)->Stress.thirdTrace() + K * pressureIncrement;
  ! }
! }


! subroutine CalcStressElastic (dt)   !!!! ONLY TEST

  ! implicit none
  ! real(fp_kind) :: SRT(3,3), RS(3,3), ident(3,3)
  ! real(fp_kind) :: strain(3,3)
  ! integer :: e,gp
  ! real(fp_kind) ,intent(in):: dt
  
  ! real(fp_kind) :: p00
  
  ! p00 = 0.
  
  ! ident = 0.0d0
  ! ident (1,1) = 1.0d0; ident (2,2) = 1.0d0; ident (3,3) = 1.0d0
  
  ! ! !!!$omp parallel do num_threads(Nproc) private (RotationRateT, Stress, SRT, RS)
  ! do e = 1, elem_count 
    ! do gp=1,elem%gausspc(e)
      ! strain = elem%str_rate(e,gp,:,:) * dt
      ! elem%sigma(e,gp,:,:) = -elem%pressure(e,gp) * ident + elem%shear_stress(e,gp,:,:)	!Fraser, eq 3.32
      ! print *, "elem ", e, ", sigma ", elem%sigma(e,gp,:,:)
    ! ! !pt%strain(i)			= dt*pt%str_rate(i + Strain;
    ! end do !gauss point
  ! end do
  ! ! !!!!$omp end parallel do    
! end subroutine CalcStressElastic


!!!!!! ONLY ELASTIC; FOR VALIDATING 
!!! ATENTION mat_G is global and mat_E not
! subroutine Calc_Elastic_Stress(dom, dt)
  ! integer :: e,gp
  ! real(fp_kind), intent (in) :: dt
  ! real(fp_kind) :: c
  ! type (dom_type), intent (in) :: dom
  
  ! c = dom%mat_E / (1.0-dom%mat_nu*dom%mat_nu)
  ! print *, "MAT C", c
  ! print *, "MAT G", mat_G
  ! do e = 1, elem_count 
    ! do gp=1,elem%gausspc(e)
      ! elem%str_tot(e,gp,:,:) = elem%str_tot(e,gp,:,:) + elem%str_inc(e,gp,:,:) * dt
      ! elem%sigma(e,gp,1,1) = c * (elem%str_tot(e,gp,1,1)+dom%mat_nu*elem%str_tot(e,gp,2,2))
      ! elem%sigma(e,gp,2,2) = c * (elem%str_tot(e,gp,2,2)+dom%mat_nu*elem%str_tot(e,gp,1,1))
      ! elem%sigma(e,gp,1,2) = 2.0* mat_G * elem%str_tot(e,gp,1,2)
      ! elem%sigma(e,gp,2,1) = elem%sigma(e,gp,1,2)
    ! end do
  ! end do  
! end subroutine

subroutine Calc_Elastic_Stress(dom, dt)
  integer :: e,gp
  real(fp_kind), intent (in) :: dt
  real(fp_kind) :: c
  type (dom_type), intent (in) :: dom
  
  !!!! PLAIN STRESS
  !c = dom%mat_E / (1.0-dom%mat_nu*dom%mat_nu)
  
  !!!! PLAIN STRAIN
  c = dom%mat_E / ((1.0+dom%mat_nu)*(1.0-2.0*dom%mat_nu))
  ! print *, "MAT C", c
  ! print *, "MAT G", mat_G
  do e = 1, elem_count 
    do gp=1,elem%gausspc(e)
      ! elem%str_inc(e,gp,:,:) = elem%str_inc(e,gp,:,:) + elem%str_inc(e,gp,:,:) * dt
      ! elem%sigma(e,gp,1,1) = elem%sigma(e,gp,1,1) + c * (elem%str_inc(e,gp,1,1)+dom%mat_nu*elem%str_inc(e,gp,2,2))
      ! elem%sigma(e,gp,2,2) = elem%sigma(e,gp,2,2) + c * (elem%str_inc(e,gp,2,2)+dom%mat_nu*elem%str_inc(e,gp,1,1))
      ! elem%sigma(e,gp,1,2) = elem%sigma(e,gp,1,2) + 2.0* mat_G * elem%str_inc(e,gp,1,2)
      ! elem%sigma(e,gp,2,1) = elem%sigma(e,gp,1,2)
      
      !!!! PLAIN STRAIN  
      elem%str_inc(e,gp,:,:) = elem%str_inc(e,gp,:,:) + elem%str_inc(e,gp,:,:) * dt
      elem%sigma(e,gp,1,1) = elem%sigma(e,gp,1,1) + c * ((1.0-dom%mat_nu)*elem%str_inc(e,gp,1,1)+dom%mat_nu*elem%str_inc(e,gp,2,2))
      elem%sigma(e,gp,2,2) = elem%sigma(e,gp,2,2) + c * ((1.0-dom%mat_nu)*elem%str_inc(e,gp,2,2)+dom%mat_nu*elem%str_inc(e,gp,1,1))
      elem%sigma(e,gp,1,2) = elem%sigma(e,gp,1,2) + (1.0-2.0*dom%mat_nu) * elem%str_inc(e,gp,1,2)
      elem%sigma(e,gp,2,1) = elem%sigma(e,gp,1,2)      
    end do
  end do 
  

end subroutine

!!!!!! IT ASSUMES PRESSURE AND STRAIN RATES ARE ALREADY CALCULATED
!!!!!! (AT t+1/2 to avoid stress at rigid rotations, see Benson 1992)
subroutine CalcStressStrain (dt) 
	use omp_lib
	
  implicit none
  !real(fp_kind) :: SRT(3,3), RS(3,3), ident(3,3
  real(fp_kind) :: SRT(3,3), RS(3,3), ident(3,3)
  integer :: e,gp
  real(fp_kind) ,intent(in):: dt
  
  real(fp_kind) :: p00, J2, sig_trial, trace, i
  
  p00 = 0.
  
  ident = 0.0d0
  ident (1,1) = 1.0d0; ident (2,2) = 1.0d0;  ident (3,3) = 1.0d0
  
  
  ! if (dim == 2) then 
    ! if (plane_mode == pl_strain) then
      
    ! end if
  ! end if
  
  ! !!!$omp parallel do num_threads(Nproc) private (RotationRateT, Stress, SRT, RS)
	!$omp parallel do num_threads(Nproc) private (e,gp,SRT,RS,J2,sig_trial) 
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
      ! if (dim ==3) then 
      trace = 0.0d0
      do i=1, 3
        trace = trace + elem%str_rate(e,gp,i,i)
      end do 
      
      elem%shear_stress(e,gp,:,:)	= dt * (2.0 * mat_G *(elem%str_rate(e,gp,:,:) - 1.0/3.0 * &
                                   (elem%str_rate(e,gp,1,1)+elem%str_rate(e,gp,2,2)+elem%str_rate(e,gp,3,3))*ident) &
                                   +SRT+RS) + elem%shear_stress(e,gp,:,:)
      
      !!! CHECK, CRASHES
      ! elem%shear_stress(e,gp,:,:)	= dt * (2.0 * mat_G *(elem%str_rate(e,gp,:,:) - 1.0/3.0 * &
                                   ! (trace)*ident) &
                                   ! +SRT+RS) + elem%shear_stress(e,gp,:,:)
      J2 = 0.5d0 * ( elem%shear_stress(e,gp,1,1) * elem%shear_stress(e,gp,1,1) + 2.0d0 * &
                      elem%shear_stress(e,gp,1,2) * elem%shear_stress(e,gp,2,1) &
                    + 2.0d0 * elem%shear_stress(e,gp,1,3) * elem%shear_stress(e,gp,3,1) &
                    + elem%shear_stress(e,gp,2,2) * elem%shear_stress(e,gp,2,2) &
                    + 2.0d0 * elem%shear_stress(e,gp,2,3) * elem%shear_stress(e,gp,3,2) &
                    + elem%shear_stress(e,gp,3,3) * elem%shear_stress(e,gp,3,3))
      
      sig_trial = sqrt(3.0d0*J2)
      !YIELD, SCALE BACK
      if (elem%sigma_y(e,gp)<sig_trial) then
        elem%shear_stress(e,gp,:,:) = elem%shear_stress(e,gp,:,:) * elem%sigma_y(e,gp) / sig_trial
        !dep=( sig_trial - sigma_y[i])/ (3.*G[i] /*+ Ep*/);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
        !pl_strain[i] += dep;	
        elem%pl_strain(e,gp) = elem%pl_strain(e,gp) + (sig_trial - elem%sigma_y(e,gp)) / (3.0d0 * mat_G) !
      endif
    
    
      ! end if
      ! print *, " shear_stress ", elem%shear_stress(e,gp,:,:)
    
      elem%sigma(e,gp,:,:) = -elem%pressure(e,gp) * ident + elem%shear_stress(e,gp,:,:)	!Fraser, eq 3.32
     
      
      ! print *, "elem ", e, ", sigma ", elem%sigma(e,gp,:,:)
      ! print *, "elem ", e, ", sigma pressure comp", -elem%pressure(e,gp)
      ! print *, "elem ", e, ", sigma shear comp", elem%shear_stress(e,gp,:,:)	
    ! !pt%strain(i)			= dt*pt%str_rate(i + Strain;
    end do !gauss point
  end do
	!$omp end parallel do

  ! !!! IF PLANE STRESS
  ! if (dim == 2) then 
    ! elem%sigma(e,gp,3,3) = 0.0
    ! elem%sigma(e,gp,3,3) = 0.0    
  ! end if 
  ! !!!!$omp end parallel do    
end subroutine CalcStressStrain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! CALC STRESS FROM ALREADY INTEGRATED STRAIN INCREMENT !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INSTEAD OF STRAIN RATE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcStress (dt) 

  implicit none
  real(fp_kind) :: SRT(3,3), RS(3,3), ident(3,3)
  integer :: e,gp
  real(fp_kind) ,intent(in):: dt
  
  real(fp_kind) :: dev(3,3)
  
  real(fp_kind) :: sig(3,3)
  
  real(fp_kind) :: p00
  
  p00 = 0.
  
  ident = 0.0d0
  ident (1,1) = 1.0d0; ident (2,2) = 1.0d0; ident (3,3) = 1.0d0
  
  ! !!!$omp parallel do num_threads(Nproc) private (RotationRateT, Stress, SRT, RS)
  do e = 1, elem_count   
    do gp=1,elem%gausspc(e)
      
      dev = deviator(elem%sigma(e,gp,:,:)) !! OBTAIN FROM DEV

      dev = dev + 2.0 * mat_G * deviator(elem%str_inc (e,gp, :,:))
      
      elem%shear_stress(e,gp,:,:) = dev

      !print *, "mat g", mat_G
      ! elem%shear_stress(e,gp,:,:)	= dt * (2.0 * mat_G *(elem%str_rate(e,gp,:,:) - 1.0/3.0 * &
                                   ! (elem%str_rate(e,gp,1,1)+elem%str_rate(e,gp,2,2)+elem%str_rate(e,gp,3,3))*ident) &
                                   ! +SRT+RS) + elem%shear_stress(e,gp,:,:)
      !print *, " shear_stress ", elem%shear_stress(e,gp,:,:)
    
      sig = dev - elem%pressure(e,gp) * ident 
      !print *, "sigma prev rot", sig
      elem%sigma(e,gp,:,:) = matmul(matmul(elem%rmat(e,gp,:,:),sig),transpose(elem%rmat(e,gp,:,:)))
      if (gp == 1) then
        print *, "sigma ", elem%sigma(e,gp,:,:)
        print *, "dev str inc ", deviator(elem%str_inc (e,gp, :,:))
        print *, "dev stress 11 22 33 12 23 31", dev(1,1),  dev(2,2),  dev(3,3), dev(1,2),  dev(2,3),  dev(3,1)
      end if 
      !!! PERFORM ROTATION
      
      !print *, "elem ", e, ", sigma ", elem%sigma(e,gp,:,:)
    ! !pt%strain(i)			= dt*pt%str_rate(i + Strain;
    end do !gauss point
  end do
  ! !!!!$omp end parallel do    
end subroutine CalcStress

subroutine calc_elem_vol ()
  implicit none
  integer :: e, gp
  real(fp_kind):: w, prev_vol, f
  
  ! P00+(Cs0*Cs0)*(Density-Density0);
  do e = 1, elem_count
    prev_vol = elem%vol(e)
    elem%vol(e) = 0.0d0
    !print *, "elem%gausspc(e)", elem%gausspc(e)
    if (elem%gausspc(e).eq.1) then
      w = 2.0**dim
    else 
      w = 1.0
    end if
    do gp=1,elem%gausspc(e)
      f = 1.0
      if (bind_dom_type .eq. 3 .and. axisymm_vol_weight .eqv. .true.) then 
        f = elem%radius(e,gp)
      endif
      !elem%vol(e) = 
      !print *, "elem e j, w", elem%detJ(e,gp), w
      elem%vol(e) = elem%vol(e) + elem%detJ(e,gp)*w*f
    end do !gp  
     !print *, "Elem ", e, "vol ",elem%vol(e)
    elem%vol_inc(e) = elem%vol(e) - prev_vol
  end do !elem

end subroutine


subroutine impose_bcv
  implicit none
  integer :: n, d
  n = 1
  do while (n .le. node_count)    
    do d=1,dim
      if (nod%is_bcv(n,d) .eqv. .true.) then
        nod%v(n,d) = nod%bcv(n,d)
        ! print *, "BCV VEL node dim", n, d 
        !print *, "nod ", n, ", ",nod%bcv(n,d), ", d", d
      end if
      
      if (nod%is_fix(n,d) .eqv. .true.) then
        ! print *, "BCV FIX node dim", n, d 
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
      if ((nod%is_fix(n,d) .eqv. .true.) .or. (nod%is_bcv(n,d).eqv. .true. )) then !!!!USE PARENTHESES
        ! print *, "BCA node dim", n, d 
        nod%a(n,d) = 0.0d0
      end if 
    end do !dim
  end do !Node    
end subroutine

subroutine CalcEquivalentStress
  use Domain
  real(fp_kind) :: J2
  integer :: e, gp
  
  do e=1,elem_count
    do gp=1,elem%gausspc(e)
    J2	= 0.5*(elem%shear_stress(e,gp,1,1)*elem%shear_stress(e,gp,1,1) + &
          2.0*elem%shear_stress(e,gp,1,2)*elem%shear_stress(e,gp,2,1) + &
          2.0*elem%shear_stress(e,gp,1,3)*elem%shear_stress(e,gp,3,1) + &
              elem%shear_stress(e,gp,2,2)*elem%shear_stress(e,gp,2,2) + &
          2.0*elem%shear_stress(e,gp,2,3)*elem%shear_stress(e,gp,3,2) + &
              elem%shear_stress(e,gp,3,3)*elem%shear_stress(e,gp,3,3));
    !print *, "J2 ", J2
    elem%sigma_eq(e,gp) = sqrt(3.0*J2)
    end do
  end do
end subroutine CalcEquivalentStress

end module
