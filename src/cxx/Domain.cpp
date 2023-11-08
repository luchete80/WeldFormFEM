#include "Domain.h"
#include "Math/Vec3D.h"
#include <iostream>ç

using namespace std;

Domain::Domain(){
  
}
void Domain::AllocateNodes(){
	
}


void Domain::AddBoxLength(Vec3D const & V, Vec3D const & L, const double &r){
	    // integer, intent(in):: tag
    // logical, intent(in) :: redint
    // !real(fp_kind), intent(in), allocatable :: V
    // real(fp_kind), dimension(1:3), intent(in)  :: V ! input
    // real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
    // real(fp_kind), dimension (1:3) :: Xp
    // integer :: i, j, k, p, ex, ey, ez, nnodz, gp
    int nodxelem;
    int nel[3];
    
    nel[0] = (int)(L(0)/(2.0*r));
    nel[1] = (int)(L(1)/(2.0*r));
    if (m_dim == 2){
      nel[2] = 0;
      nodxelem = 4;
    } else {
      nel[2] = (int)(L(2)/(2.0*r));
      nodxelem = 8;
    }
    
		Vec3D Xp;
    //Xp(2) = V(2) ;
    

    // write (*,*) "Creating Mesh ...", "Elements ", nel(1), ", ",nel(2)
    
    // call AllocateNodes((nel(1) +1)* (nel(2)+1) * (nel(3)+1))
    // print *, "Element count in XYZ: ", nel(:)
    // write (*,*) "Box Node count ", node_count
    
    
    // write (*,*) "xp ", Xp(:)    
    
    if (m_dim == 2) {
    cout << "Box Particle Count is " << node_count <<endl;
    // p = 1
    // !do while (Xp(3) <= (V(3)+Lz))
      // j = 1;         Xp(2) = V(2)
      // do while (j <= (nel(2) +1))
        // i = 1
        // Xp(1) = V(1)
        // do while (i <= (nel(1) +1))
          // nod%x(p,:) = Xp(:)
          // print *,"node ",p , "X: ",Xp(:)
          // p = p + 1
          // Xp(1) = Xp(1) + 2 * r
          // i = i +1
        // end do
        // Xp(2) = Xp(2) + 2 * r
        // j = j +1
      // end do 
      // Xp(3) = Xp(3) + 2 * r
    // !end do
    
    } else {
      // p = 1
      // k = 1; Xp(3) = V(3)
      // do while (k <= (nel(3) +1))
        // j = 1;         Xp(2) = V(2)
        // do while (j <= (nel(2) +1))
          // i = 1
          // Xp(1) = V(1)
          // do while (i <= (nel(1) +1))
            // nod%x(p,:) = Xp(:)
            // print *,"node ",p , "X: ",Xp(:)
            // p = p + 1
            // Xp(1) = Xp(1) + 2 * r
            // i = i +1
          // end do
          // Xp(2) = Xp(2) + 2 * r
          // j = j +1
        // end do 
        // Xp(3) = Xp(3) + 2 * r
        // k = k + 1
      // end do    
		}//if dim
    
    // !! ALLOCATE ELEMENTS
    // !! DIMENSION = 2
    // gp = 1
    // if (dim .eq. 2) then
      // if (redint .eqv. .False.) then
        // gp = 4
      // end if 
      // call AllocateElements(nel(1) * nel(2),gp) !!!!REDUCED INTEGRATION
    // else 
      // if (redint .eqv. .False.) then
        // gp = 8
      // end if 
      // call AllocateElements(nel(1) * nel(2)*nel(3),gp) 
    // end if
    
    // if (dim .eq. 2) then
      // ey = 0
      // i = 1
      // do while ( ey < nel(2))
          // ex = 0
          // do while (ex < nel(1)) 
              // elem%elnod(i,:)=[(nel(1)+1)*ey + ex+1,(nel(1)+1)*ey + ex+2,(nel(1)+1)*(ey+1)+ex+2,(nel(1)+1)*(ey+1)+ex+1]         
              // print *, "Element ", i , "Elnod", elem%elnod(i,:) 
              // i=i+1
            // ex = ex + 1
          // end do
        // ey = ey + 1
      // end do  
    // else 
      // ez = 0
      // i = 1
      // nnodz = (nel(1)+1)*(nel(2)+1)
      // print *, "Element Nodes at z ", nnodz
      // do while ( ez < nel(3))
        // ey = 0    
        // do while ( ey < nel(2))
            // ex = 0
            // do while (ex < nel(1)) 
                // !elem%elnod(i,:)=[(nel(1)+1)*(ey+1)+ex+2,(nel(1)+1)*(ey+1)+ex+1,(nel(1)+1)*ey + ex+1,(nel(1)+1)*ey + ex+2]       
                // elem%elnod(i,:) = [ nnodz*ez + (nel(1)+1)*ey + ex+1,nnodz*ez + (nel(1)+1)*ey + ex+2, &
                                    // nnodz*ez + (nel(1)+1)*(ey+1)+ex+2,nnodz*ez + (nel(1)+1)*(ey+1)+ex+1, &
                                    // nnodz*(ez + 1) + (nel(1)+1)*ey + ex+1,nnodz*(ez + 1) + (nel(1)+1)*ey + ex+2, &
                                    // nnodz*(ez + 1) + (nel(1)+1)*(ey+1)+ex+2,nnodz*(ez + 1)+ (nel(1)+1)*(ey+1)+ex+1]
                // print *, "Element ", i , "Elnod", elem%elnod(i,:) 
                // i=i+1
              // ex = ex + 1
            // end do
          // ey = ey + 1
        // end do 
        // ez = ez + 1
      // end do !el z
    // end if !dim
    
    // call AllocateDomain()
    // i = 1
    // do while ( i <= node_count)
      // nod%is_bcv(i,:) = .false.
      // i = i + 1
    // end do
  
    // ! nod%m(:)   = Density * Lx * Ly * Lz / node_count
    // ! nod%rho(:)   = Density
    // elem%rho_0(:,:) = Density
    // !print *, "Particle mass ", nod%m(2)
    
    // !nod%id(:) = tag
    
    // fext_glob (:,:) = 0.0d0
    
    // elem%e_length(:) = Lx !TODO: CHANGE!
    
    // tot_mass = Density * Lx * Ly * Lz
    // if (dim == 2) then !!!assuming plain strain
      // tot_mass = tot_mass / Lz
    // end if
    // print *, "Total Mass: ", tot_mass
    
    // call SearchNodelem
}