#include "Domain.cuh"

// void Domain:: Calc_elem_strains (){
  // // implicit none
  // // integer :: e, i,j,k, gp, d, n
  // // real(fp_kind), dimension(dim,nodxelem) ::temp
  // // real(fp_kind) :: f
  // // real(fp_kind) :: test(1,6) !ifwanted to test in tensor form
  
  // // elem%str_rate = 0.0d0
  // // elem%rot_rate = 0.0d0
  
  // for (int e=0;elem_count;e++){
    // // do gp = 1, elem%gausspc(e)
      // // !Is only linear matrix?    
      // // !!!TODO: CHANGE FROM MATRIX OPERATION TO SIMPLE OPERATION
      // // f = 1.0d0/elem%detJ(e,gp)
      // // temp = elem%dHxy_detJ(e,gp,:,:) * f!!!!TODO: MODIFY BY MULTIPLYING
      // // elem%strain(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%uele (e,:,:)) 
      // // !print *, "standard stran rate calc (matricial) "
      // // ! !!!! DEFAULT (TODO: CHECK IF IS SLOW)
      // // !test = f* matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:))  ! (6x24)(24x1)
      // // !print *, "e11 e22 e33 2e12 2e23 2e31", test

      // // do n=1, nodxelem  
        // // do d=1, dim
          // // !print *, "node dim dHxy vele", n,d,temp(d,n) , elem%vele (e,dim*(n-1)+d,1) 
          // // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // // elem%rot_rate(e,gp, d,d) = 0.0d0
        // // end do
        // // !!!! TO AVOID ALL MATMULT
        // // elem%str_rate(e,gp, 1,2) = elem%str_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) &!!!!dvx/dy
                                   // // + temp(1,n) * elem%vele (e,dim*(n-1)+2,1)
        // // elem%rot_rate(e,gp, 1,2) = elem%rot_rate(e,gp, 1,2) + temp(2,n)* elem%vele (e,dim*(n-1)+1,1) & !!!!dvx/dx
                                   // // - temp(1,n) * elem%vele (e,dim*(n-1)+2,1)                           !!!!
        // // if (dim == 3) then
          // // elem%str_rate(e,gp, 2,3) = elem%str_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &!!!d/dz*vy     
                                     // // + temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // // elem%str_rate(e,gp, 1,3) = elem%str_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // // + temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz     
          // // elem%rot_rate(e,gp, 2,3) = elem%rot_rate(e,gp, 2,3) + temp(3,n)* elem%vele (e,dim*(n-1)+2,1) &
                                     // // - temp(2,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dy*vz
          // // elem%rot_rate(e,gp, 1,3) = elem%rot_rate(e,gp, 1,3) + temp(3,n)* elem%vele (e,dim*(n-1)+1,1) & !!!d/dz*vx     
                                     // // - temp(1,n) * elem%vele (e,dim*(n-1)+3,1)    !!!d/dx*vz    
        // // end if     
      // // end do !Nod x elem
      // // elem%str_rate(e,gp, 1,2) = 0.5 * elem%str_rate(e,gp, 1,2); 
      // // elem%rot_rate(e,gp, 1,2) = 0.5 * elem%rot_rate(e,gp, 1,2)      

      // // elem%str_rate(e,gp, 2,1) =     elem%str_rate(e,gp, 1,2)
      // // elem%rot_rate(e,gp, 2,1) =    -elem%rot_rate(e,gp, 1,2)
      // // if (dim .eq. 3) then
        // // elem%str_rate(e,gp, 1,3) = 0.5 * elem%str_rate(e,gp, 1,3); elem%str_rate(e,gp, 2,3) = 0.5 * elem%str_rate(e,gp, 2,3)
        // // elem%rot_rate(e,gp, 1,3) = 0.5 * elem%rot_rate(e,gp, 1,3); elem%rot_rate(e,gp, 2,3) = 0.5 * elem%rot_rate(e,gp, 2,3)
        
        // // elem%str_rate(e,gp, 3,2) =     elem%str_rate(e,gp, 2,3)
        // // elem%str_rate(e,gp, 3,1) =     elem%str_rate(e,gp, 1,3)

        // // elem%rot_rate(e,gp, 3,2) =     -elem%rot_rate(e,gp, 2,3)
        // // elem%rot_rate(e,gp, 3,1) =     -elem%rot_rate(e,gp, 1,3)
      // // end if

      // // !elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
      // // !print *, "simlpified strain rate "
      // // !print *, "strain rate ", elem%str_rate(e,gp,:,:)
      // // !print *, "rot    rate ", elem%rot_rate(e,gp,:,:)
    // // end do !gp
  // }
// }


// void Domain:: Cal_elem_forces ()
  // implicit none
  // integer :: e, i,j,k, gp,n, d
  // real(fp_kind), dimension(dim*nodxelem,1) ::f
  // real(fp_kind) :: w
  // !TESTING
  // real (fp_kind) :: sigma_test(6,1) !ORDERED
  // real(fp_kind) :: test(24,1) !ifwanted to test in tensor form
  // elem%f_int = 0.0d0
  // w = 1.0d0 !!! Full integration

  // do e=1, elem_count
    // if (elem%gausspc(e) .eq. 1) then
      // w = 2.0d0**dim
    // end if
    // do gp = 1, elem%gausspc(e)
      // !print *, "elem%dHxy_detJ(e,gp,1", elem%dHxy_detJ(e,gp,1,:)
      // !print *, "elem%dHxy_detJ(e,gp,2", elem%dHxy_detJ(e,gp,2,:)
      // sigma_test (:,1)=[elem%sigma (e,gp, 1,1),elem%sigma (e,gp, 2,2),elem%sigma (e,gp, 3,3),&
                        // elem%sigma (e,gp, 1,2),elem%sigma (e,gp, 2,3),elem%sigma (e,gp, 3,1)]
      // test = w*matmul(transpose(elem%bl(e,gp,:,:)),sigma_test)  ! (24x6)(6x1)
      // !print *, "test force", test
      
      // print *, "dHdxy, 1", elem%dHxy_detJ(e,gp,1,:)
      // print *, "dHdxy, 2", elem%dHxy_detJ(e,gp,2,:)
      // !print *, "dHdxy, 3", elem%dHxy_detJ(e,gp,1,:)
      
      // do n=1, nodxelem
        // do d=1, dim
          // elem%f_int(e,n,d) = elem%f_int(e,n,d) + elem%dHxy_detJ(e,gp,d,n) * elem%sigma (e,gp, d,d)
        // end do
        // if (dim .eq. 2) then  !!!!! TODO: CHANGE WITH BENSON 1992 - EQ 2.4.2.11 FOR SIMPLICITY
          // elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) 
          // elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2)
        // else 
          // elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) + &
                                                  // elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 1,3)
          // elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) + &
                                                  // elem%dHxy_detJ(e,gp,3,n) * elem%sigma (e,gp, 2,3)
          // elem%f_int(e,n,3) = elem%f_int(e,n,3) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 2,3) + &
                                                  // elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,3)
        // end if
        // print *, "Element force Node ", n, "F  ", elem%f_int(e,n,:) * w 
      // end do! nod x elem
      // !print *, "test ", w * elem%dHxy_detJ(e,gp,3,8)  * elem%sigma (e,gp, 3,3)
      // !print *, "dHxy ", elem%dHxy_detJ(e,gp,3,8), "w ", w
      // !print *, "s33 ", elem%sigma (e,gp, 3,3)
    // end do !gp
    // elem%f_int(e,:,:) = elem%f_int(e,:,:) * w
  // end do!elem
// end subroutine