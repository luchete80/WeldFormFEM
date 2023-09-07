
// // !!!! IN ORDER TO CALCULATE IT ONLY ONCE
// // !!!!!!!!!!!! IN RS PLANE ---->>>>
// // !!!!!!!! 4------3
// // !!!!!!!! |      |
// // !!!!!!!! 1 ---- 2
// // !!!! CALCULATE JACOBIAN AND DETERMINANT
#include "ElementData.cuh"
#include "NodeData.cuh"
#include "Math.cuh"

__device__ inline void calculate_element_Jacobian(NodeData *nod, ElementData *elem, const int &dim, const int &gp) {
  
  int e = threadIdx.x + blockDim.x*blockIdx.x;
// subroutine calculate_element_Jacobian ()
  // integer :: e
  // ! !rg=gauss[ig]
  // ! !sg=gauss[jg]
  // real(fp_kind), dimension(dim,nodxelem) :: dHrs !!! USED ONLY FOR SEVERAL GAUSS POINTS
  // real(fp_kind), dimension(nodxelem,dim) :: x2
  // real(fp_kind), dimension(dim,dim) :: test
  // real(fp_kind), dimension(dim, dim*nodxelem) :: temph
  
  // integer :: i,j,k, gp
  // real(fp_kind):: r   !!! USED ONLY FOR SEVERAL GAUSS POINTS
  // real(fp_kind), dimension(8,3):: gpc !!! gauss point coordinates, r,s,t
  
  // gp = 1
  // do e=1, elem_count
// ! #ifdef _PRINT_DEBUG_  
    // ! print *, "el ", e 
// ! #endif    
    // do i=1,nodxelem
        // !print *, "elnod " , elem%elnod(e,i)
        // x2(i,:)=nod%x(elem%elnod(e,i),:)
    // end do
    
    if (elem->gausspc[e] == 1) {      
    
      // if (dim == 2) then 
        // !dHdrs [-1,1,1,-1;  -1.-1,1,1] x X2
        // !! J = [
        // !! dx/dr dy/dr
        // !! dx/ds dy/dx ]
        // !!! THIS IS TO AVOID MATMUL
        // ! print *, "nodes X ", x2(:,1)
        // ! print *, "nodes Y ", x2(:,2)
                
        // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)
        // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)
        // elem%jacob(e,gp,:,:) = 0.25*elem%jacob(e,gp,:,:)
        // else !!!DIM 3
          // !!!!! SETTING LIKE THIS AVOID MATMUL
          // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)-x2(5,:)+x2(6,:)+x2(7,:)-x2(8,:)
          // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)-x2(5,:)-x2(6,:)+x2(7,:)+x2(8,:)
          // elem%jacob(e,gp,3,:) = -x2(1,:)-x2(2,:)-x2(3,:)-x2(4,:)+x2(5,:)+x2(6,:)+x2(7,:)+x2(8,:)
          // !elem%jacob(e,gp,2,:) = [-x2(1,2),-x2(2,2), x2(3,2), x2(4,2),-x2(5,2),-x2(6,2), x2(7,2), x2(8,2)]
          // !elem%jacob(e,gp,3,:) = [-x2(1,3),-x2(2,3), x2(3,3), x2(4,3),-x2(5,3),-x2(6,3), x2(7,3), x2(8,3)]
          // ! dHrs(1,:)=[-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0] AND THIS IS dHrs*x2
          // ! dHrs(2,:)=[-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0]       
          // ! dHrs(3,:)=[-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0]  
          // ! elem%jacob(e,gp,1,:) = matmul(dHrs,x2)
          // elem%jacob(e,gp,:,:) = 0.125*elem%jacob(e,gp,:,:)
      // end if  !!!!DIM
      // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
    } else { //!!!!! GP > 1
      double r = 1.0/sqrt(3.0);
      double gpc[8][3] = {  {-r,-r,-r},   {r,-r,-r},      {-r, r,-r}, {r, r,-r}, //!These are the 4 points for 2D full elem
                          {-r,-r, r},   {r,-r, r},      {-r, r, r}, {r, r, r}};
    
      if (dim == 3) {
        // do gp = 1,8
          // dHrs(1,:)=[-1.0*(1-gpc(gp,2))*(1.0-gpc(gp,3)),     (1-gpc(gp,2))*(1.0-gpc(gp,3))&
                    // ,     (1+gpc(gp,2))*(1.0-gpc(gp,3)),-1.0*(1+gpc(gp,2))*(1.0-gpc(gp,3))&
                    // ,-1.0*(1-gpc(gp,2))*(1.0+gpc(gp,3)),     (1-gpc(gp,2))*(1.0+gpc(gp,3))&
                    // ,     (1+gpc(gp,2))*(1.0+gpc(gp,3)),-1.0*(1+gpc(gp,2))*(1.0+gpc(gp,3))]
          // dHrs(2,:)=[-1.0*(1-gpc(gp,1))*(1.0-gpc(gp,3)),-1.0*(1+gpc(gp,1))*(1.0-gpc(gp,3))&
                         // ,(1+gpc(gp,1))*(1.0-gpc(gp,3)),     (1-gpc(gp,1))*(1.0-gpc(gp,3))&
                    // ,-1.0*(1-gpc(gp,1))*(1.0+gpc(gp,3)),-1.0*(1+gpc(gp,1))*(1.0+gpc(gp,3))&
                         // ,(1+gpc(gp,1))*(1.0+gpc(gp,3)),     (1-gpc(gp,1))*(1.0+gpc(gp,3))]
          // dHrs(3,:)=[-1.0*(1-gpc(gp,1))*(1.0-gpc(gp,2)),-1.0*(1+gpc(gp,1))*(1.0-gpc(gp,2))&
                    // ,-1.0*(1+gpc(gp,1))*(1.0+gpc(gp,2)),-1.0*(1-gpc(gp,1))*(1.0+gpc(gp,2))&
                    // ,     (1-gpc(gp,1))*(1.0-gpc(gp,2)),     (1+gpc(gp,1))*(1.0-gpc(gp,2))&
                    // ,     (1+gpc(gp,1))*(1.0+gpc(gp,2)),     (1-gpc(gp,1))*(1.0+gpc(gp,2))]                     
          
          // elem%dHrs(e,gp,:,:) =  dHrs(:,:)         
          // !dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]         
          // !dHrs(3,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))] 
          // !print *, "dhrs", dHrs 
          // !print *, "x2", x2 
          //elem%jacob(e,gp,:,:) = 0.125*matmul(dHrs,x2)
          //elem->jacob(e*gp)
// ! #if defined _PRINT_DEBUG_
          // ! print *, "jacob ", elem%jacob(e,gp,:,:)
// ! #endif          
          // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
          // !print *, "detJ ", elem%detJ(e,gp)
        // end do !gp
      } else { //!dim =2
        // do gp = 1,4
          // dHrs(1,:)=[-1.0*(1-gpc(gp,2)),     (1-gpc(gp,2))&
                    // ,     (1+gpc(gp,2)),-1.0*(1+gpc(gp,2))]
          // dHrs(2,:)=[-1.0*(1-gpc(gp,1)),-1.0*(1+gpc(gp,1))&
                         // ,(1+gpc(gp,1)),     (1-gpc(gp,1))]                
          
          // elem%dHrs(e,gp,:,:) =  dHrs(:,:)         
          // !dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]         
          // !dHrs(3,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))] 
          // !print *, "dhrs", dHrs 
          // !print *, "x2", x2 
          // elem%jacob(e,gp,:,:) = 0.25*matmul(dHrs,x2)
// ! #if defined _PRINT_DEBUG_
          // !print *, "jacob ", elem%jacob(e,gp,:,:)
// ! #endif          
          // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
          // !print *, "detJ ", elem%detJ(e,gp)
        // end do !gp      
        
      }// end if !dim
    } // end if !!gp ==1
    
// ! #if defined _PRINT_DEBUG_
    // !print *, "jacob ", elem%jacob(e,gp,:,:)
// ! #endif    
  // end do !element

}

// // // !!!!! ASUMES CALCULATED DETJ AND LOCAL EDRIV MATRICES (dHrs)
__device__ inline void calculate_element_derivMat(NodeData *nod, ElementData *elem, const int &dim, const int &gp) {
  
  int e = threadIdx.x + blockDim.x*blockIdx.x;

  // integer :: e,d
  // ! !rg=gauss[ig]
  // ! !sg=gauss[jg]
  // real(fp_kind), dimension(dim,nodxelem) :: dHrs
  // real(fp_kind), dimension(dim,dim) :: test, invJ
  
  // integer :: i,j,k, gp
  // real(fp_kind), dimension(2) :: r, s

  // !! Update x2 vector (this is useful for strain and stress things)
  

    // gp = 1
    if (elem->gausspc[e] == 1) {
      
      // !!!!!!invJ = adj(elem%jacob(e,gp,:,:))/elem%detJ(e,gp) !!!! ALREADY CALCULATED   
      // invJ = adj(elem%jacob(e,gp,:,:)) !!! IN FACT IS invJ x detJ
      if (dim == 2) {       
          // !!!!! J-1 = dr/dx
          // !!!! dHxy = J-1 x dHrs = [ ] x 0.25[-1 1 -1 1]
          // !!!!                               [-1 -1 1 1]
          // elem%dHxy_detJ(e,gp,:,1) = -invJ(:,1)-invJ(:,2) !For each 3 rows of inv J and dHdxy
          // elem%dHxy_detJ(e,gp,:,2) =  invJ(:,1)-invJ(:,2)
          // elem%dHxy_detJ(e,gp,:,3) =  invJ(:,1)+invJ(:,2)
          // elem%dHxy_detJ(e,gp,:,4) = -invJ(:,1)+invJ(:,2)     
          
          // elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.25d0
      } else { //!!!DIM 3
          
          // !!!! REMAINS ELEM_VOLUME
          // !!!!! J-1 = dr/dx
          // !!!! dHxy = J-1 x dHrs 
          // elem%dHxy_detJ(e,gp,:,1) = -invJ(:,1)-invJ(:,2)-invJ(:,3) !For each 3 rows of inv J and dHdxy
          // elem%dHxy_detJ(e,gp,:,2) =  invJ(:,1)-invJ(:,2)-invJ(:,3)
          // elem%dHxy_detJ(e,gp,:,3) =  invJ(:,1)+invJ(:,2)-invJ(:,3)
          // elem%dHxy_detJ(e,gp,:,4) = -invJ(:,1)+invJ(:,2)-invJ(:,3)
          // elem%dHxy_detJ(e,gp,:,5) = -invJ(:,1)-invJ(:,2)+invJ(:,3)
          // elem%dHxy_detJ(e,gp,:,6) =  invJ(:,1)-invJ(:,2)+invJ(:,3)
          // elem%dHxy_detJ(e,gp,:,7) =  invJ(:,1)+invJ(:,2)+invJ(:,3)
          // elem%dHxy_detJ(e,gp,:,8) = -invJ(:,1)+invJ(:,2)+invJ(:,3)

          // elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.125d0          
          // ! print *, "1,1", elem%dHxy(e,gp,1,1), "inv 1,1 1,2 1,3", invJ(1,1), invJ(1,2), invJ(1,3)
          // ! print *, "1,1", elem%dHxy(e,gp,2,1), "inv 1,1 1,2 1,3", invJ(2,1), invJ(2,2), invJ(2,3)
          // do k=1,nodxelem !!! TABLE 6.6 page 556 bathe, 24 x 6
            // do d=1,dim
              // elem%bl(e,gp,d,dim*(k-1)+d ) = elem%dHxy_detJ(e,gp,d,k) 
            // end do
            // elem%bl(e,gp,4,dim*(k-1)+1) = elem%dHxy_detJ(e,gp,2,k)
            // elem%bl(e,gp,4,dim*(k-1)+2) = elem%dHxy_detJ(e,gp,1,k)
            
            // elem%bl(e,gp,5,dim*(k-1)+2) = elem%dHxy_detJ(e,gp,3,k)
            // elem%bl(e,gp,5,dim*(k-1)+3) = elem%dHxy_detJ(e,gp,2,k)
            
            // elem%bl(e,gp,6,dim*(k-1)+1) = elem%dHxy_detJ(e,gp,3,k)
            // elem%bl(e,gp,6,dim*(k-1)+3) = elem%dHxy_detJ(e,gp,1,k)         
          // end do
      }// end if  !!!!DIM
    } else { //!!!!! GP > 1
      // ! if (dim .eq. 2) then 

      if (dim == 3){ //then  !!!! dim 3
        // do gp = 1,8
          // invJ = adj(elem%jacob(e,gp,:,:))!!!!/elem%detJ(e,gp) !!!! ALREADY CALCULATED    
          // !print *, "detJ", elem%detJ(e,gp)
          // !print *, "invJ", invJ
          // elem%dHxy_detJ(e,gp,:,:) = 0.125d0 * matmul(invJ,elem%dHrs(e,gp,:,:))
          // !print *, "dHxy", elem%dHxy_detJ(e,gp,:,:)
        // end do !!!!gp
      // else !dim =2 
        // do gp = 1,4
          // invJ = adj(elem%jacob(e,gp,:,:))!!!!/elem%detJ(e,gp) !!!! ALREADY CALCULATED    
          // !print *, "detJ", elem%detJ(e,gp)
          // !print *, "invJ", invJ
          // elem%dHxy_detJ(e,gp,:,:) = 0.25d0 * matmul(invJ,elem%dHrs(e,gp,:,:))
          // !print *, "dHxy", elem%dHxy_detJ(e,gp,:,:)
        // end do !!!!gp      
    }// end if!!!di 3
  }// end if
}