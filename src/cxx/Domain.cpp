#include "Domain.h"
#include "Math/Vec3D.h"
#include <iostream>
#include "Node.h"
#include "El4N2DPE.h"

using namespace std;

Domain::Domain(){
  
}
void Domain::AllocateNodes(const int &nc){
	node_count = nc;
}


void Domain::AddBoxLength(Vec3D const & V, Vec3D const & L, const double &r){
	    // integer, intent(in):: tag
    // logical, intent(in) :: redint
    // !real(fp_kind), intent(in), allocatable :: V
    // real(fp_kind), dimension(1:3), intent(in)  :: V ! input
    // real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
    Vec3D Xp;
    int p, nnodz;
    int nodxelem;
    int nel[3];
    
    nel[0] = (int)(L(0)/(2.0*r));
    nel[1] = (int)(L(1)/(2.0*r));
    cout << "Nel x: "<<nel[0]<<", y "<<nel[1]<<endl;
    if (m_dim == 2){
      nel[2] = 0;
      nodxelem = 4;
    } else {
      nel[2] = (int)(L(2)/(2.0*r));
      nodxelem = 8;
    }
    

    Xp(2) = V(2) ;
    

    // write (*,*) "Creating Mesh ...", "Elements ", nel(1), ", ",nel(2)
    
    AllocateNodes((nel[0] +1) * (nel[1]+1) * (nel[2]+1));
    // print *, "Element count in XYZ: ", nel(:)
    // write (*,*) "Box Node count ", node_count
    
    
    // write (*,*) "xp ", Xp(:)    
    
    if (m_dim == 2) {
    cout << "Box Particle Count is " << node_count <<endl;
    p = 1;
      Xp(1) = V(1);
      for (int j = 0; j < (nel[1] +1);j++){
        Xp(0) = V(0);
        for (int i = 0; i < (nel[0] +1);i++){
					m_node.push_back(new Node(Xp));
          //nod%x(p,:) = Xp(:);
          cout << "node " << p <<"X: "<<Xp<<endl;
          p++;
          Xp(0) = Xp(0) + 2 * r;
        }
        Xp(1) = Xp(1) + 2 * r;
      }// 
      Xp(2) = Xp(2) + 2 * r;

    cout <<"m_node size"<<m_node.size()<<endl;
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
    int gp = 1;
    if (m_dim == 2) {
      // if (redint .eqv. .False.) then
        // gp = 4
      // end if 
      //call AllocateElements(nel(1) * nel(2),gp) !!!!REDUCED INTEGRATION
    } else {
      // if (redint .eqv. .False.) then
        // gp = 8
      // end if 
      // call AllocateElements(nel(1) * nel(2)*nel(3),gp) 
    }
    
		int ex, ey, ez;
		std::vector<Node*> n;
		
    if (m_dim == 2) {
			n.resize(4);
      for (int ey = 0; ey < nel[1];ey++){

           for (int ex = 0; ex < nel[0];ex++){
        int iv[4];
        iv[0] = (nel[0]+1)*ey + ex;        iv[1] = (nel[0]+1)*ey + ex+1;
        iv[2] = (nel[0]+1)*(ey+1) + ex+1;        iv[3] = (nel[0]+1)*(ey+1) + ex;
        // cout << i[]
						n[0]= m_node[iv[0]];
						n[1]= m_node[(nel[0]+1)*ey + ex+1];
						n[2]= m_node[(nel[0]+1)*(ey+1)+ex+1];
						n[3]= m_node[(nel[0]+1)*(ey+1)+ex];
            cout << "Nel x : "<<nel[0]<<endl;
           cout << "nodes "<<endl;
           for (int i=0;i<4;i++)cout << iv[i]<<", ";
						 m_element.push_back(new El4N2DPE(n));
																							// m_node[(nel[0]+1)*ey + ex+1],
																							// m_node[(nel[0]+1)*(ey+1)+ex+1],
																							// m_node[(nel[0]+1)*(ey+1)+ex]
																							// );
              //elem%elnod(i,:)=[(nel(1)+1)*ey + ex+1,(nel(1)+1)*ey + ex+2,(nel(1)+1)*(ey+1)+ex+2,(nel(1)+1)*(ey+1)+ex+1]         
              //print *, "Element ", i , "Elnod", elem%elnod(i,:) 
					 }
      } 
    } else {
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
		}//if dim 
    
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