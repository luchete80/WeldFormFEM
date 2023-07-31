#include "Domain.cuh"
#include <iostream>

using namespace std;

void Domain::AddBoxLength(const int &tag, double3 V, double Lx, double Ly, double Lz, double &r){		
    // integer, intent(in):: tag
    // logical, intent(in) :: redint
    // !real(fp_kind), intent(in), allocatable :: V
    // real(fp_kind), dimension(1:3), intent(in)  :: V ! input
    // real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
    double3 Xp;
    int i, j, k, p, ex, ey, ez, nnodz, gp;
      
    int nel [3];
    
    
    nel[0] = (int)(Lx/(2.0*r)); 
    nel[1] = (int)(Ly/(2.0*r)); 
    if (dim == 2) {
      nel[2] = 0;
      nodxelem = 4;
    } else {
      nel [2] = (int)(Lz/(2.0*r)); 
      nodxelem = 8;
    }
    
    // Xp(3) = V(3) 
    

    cout << "Creating Mesh ..."<< "Elements "<< nel[0]<< ", "<<nel[1]<<endl;
    int nc = (nel[0] +1)* (nel[1]+1) * (nel[2]+1);
    cout << "Nodes: "<<nc<<endl;
    
    AllocateNodes(nc);
    
    double3 *x_h =  new double3 [nc];
    // print *, "Element count in XYZ: ", nel(:)
    // write (*,*) "Box Node count ", node_count
    
    
    // write (*,*) "xp ", Xp(:)    
    
    //Allocate nodes
    if (dim ==2) {
    // !write(*,*) "Box Particle Count is ", node_count
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
    // end do
    
    } else {
      p = 0;
      Xp.z = V.z;
      for (k=1;k <= (nel[2] +1);k++){
        Xp.y  = V.y;
        for (j=1;j <= (nel[1] +1);j++){
          Xp.x = V.x;
          for (i=1;i <= (nel[0] +1);i++){
            x_h[p] = Xp;
            //print *,"node ",p , "X: ",Xp[:]
            cout << x_h[p].x<<", "<<x_h[p].y<<", "<<x_h[p].z<<endl;
            p = p + 1;
            Xp.x += 2.0 * r;
          }
          Xp.y += 2.0 * r;
        }
        Xp.z += 2.0 * r;
      }//k  
    }//dim 3
    
    cudaMemcpy(this->x, x_h, nc * sizeof(double3), cudaMemcpyHostToDevice);	
    
    int ne;
    // !! ALLOCATE ELEMENTS
    // !! DIMENSION = 2
    gp = 1;
    if (dim == 2) {
      
      if (!redint) gp = 4;
      ne = nel[0] * nel[1];
      nodxelem = 4;
    }else {
      if (!redint) gp = 8;
      ne = nel[0] * nel[1]*nel[2];
      nodxelem = 8;
    }
    AllocateElements(ne,gp);
    
    unsigned long *elnod_x =  new unsigned long [ne * nodxelem];
    if (dim == 2) {
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
    } else {
      ez = 0;
      i = 0;
      nnodz = (nel[0]+1)*(nel[1]+1);
      cout << "Element Nodes at z " << nnodz;
      for (ez=0; ez < nel[2];ez++){
        ey = 0   ;
        for (ey=0;ey < nel[1];ey++){
            for (ex=0;ex < nel[0];ex++){ 
                //!elem%elnod(i,:)=[(nel(1)+1)*(ey+1)+ex+2,(nel(1)+1)*(ey+1)+ex+1,(nel(1)+1)*ey + ex+1,(nel(1)+1)*ey + ex+2]       
                
                elnod_x[i*nodxelem  ] = nnodz*ez + (nel[0]+1)*ey + ex;
                elnod_x[i*nodxelem+1] = nnodz*ez + (nel[0]+1)*ey + ex+1;
                elnod_x[i*nodxelem+2] = nnodz*ez + (nel[0]+1)*(ey+1)+ex+1;
                elnod_x[i*nodxelem+3] = nnodz*ez + (nel[0]+1)*(ey+1)+ex;
                elnod_x[i*nodxelem+4] = nnodz*(ez + 1) + (nel[0]+1)*ey + ex;
                elnod_x[i*nodxelem+5] = nnodz*(ez + 1) + (nel[0]+1)*ey + ex+1;
                elnod_x[i*nodxelem+6] = nnodz*(ez + 1) + (nel[0]+1)*(ey+1)+ex+1;
                elnod_x[i*nodxelem+7] = nnodz*(ez + 1) + (nel[0]+1)*(ey+1)+ex;
                cout << "Element "<< i << "Elnod"<<endl;
                for (int d=0;d<nodxelem;d++) cout << elnod_x[i*nodxelem+d]<<", ";
                cout <<endl;
                i++;
            }//ex
        } //ey
      }//ez
    }
    //Already allocated
    cudaMemcpy(this->elnod, elnod, elem_count * nodxelem * sizeof (unsigned long ), cudaMemcpyHostToDevice);	
    
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
    
    // tot_mass = Density * Lx * Ly * Lz
    // print *, "Total Mass: ", tot_mass
    
    // call SearchNodelem
  // end subroutine AddBoxLength
  
 }
 
 void Domain::AllocateNodes(const int &nc){
   node_count = nc;
  cudaMalloc((void **)&x, node_count * sizeof (double3));
  cudaMalloc((void **)&v, node_count * sizeof (double3));
  cudaMalloc((void **)&a, node_count * sizeof (double3));
  cudaMalloc((void **)&u, node_count * sizeof (double3));
 }
 
 //Assumed nodxelem is already set
  void Domain::AllocateElements(const int &ne, const int &gp){
  elem_count = ne;
  cudaMalloc((void **)&elnod, elem_count * nodxelem * sizeof (unsigned long ));
  // cudaMalloc((void **)&v, node_count * sizeof (double3));
  // cudaMalloc((void **)&a, node_count * sizeof (double3));
  // cudaMalloc((void **)&u, node_count * sizeof (double3));
 }
 