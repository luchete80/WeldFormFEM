#include "Domain_d.h"
#include <iostream>
#include <vector>

#include "Matrix_temp.h"
#include <sstream>
#include <fstream> 
#include <iostream>

#ifdef CUDA_BUILD

#include "vector_math.h"
#include "tensor.cuh"

#else
// #include "Matrix.h"
#endif 

#if CUDA_BUILD
#include "tensor.cuh"
#else
//#include "Matrix.h"
#endif

#include "tensor3.C"

using namespace std;
using namespace LS_Dyna;

namespace MetFEM {

void Domain_d::SetDimension(const int &node_count, const int &elem_count){
  
  m_node_count = node_count;
  m_elem_count = elem_count;
  
  // NODAL VARIABLES
  
  malloc_t (x,      double,node_count*3);
  malloc_t (v,      double,node_count*3);
  malloc_t (a,      double,node_count*3);
  malloc_t (u,      double,node_count*3);
  malloc_t (u_dt,   double,node_count*3);
  
  malloc_t (prev_a, double,node_count*3);  
	//cudaMalloc((void **)&m_f, node_count * sizeof (double) * 3);
  malloc_t (m_fi,double,node_count*3); //Internal forces
  malloc_t (m_fe,double,node_count*3);
  
  malloc_t (m_mdiag, double,node_count);
  malloc_t (m_mglob, double,node_count*node_count); //TODO: MAKE SPARSE. DEALLOCATED AFER DIAG CALCULATION
	
  /// MATRICES ///
  /// dHxy_detJ: DIM X NODXELEM
  /// TO AVOID EXCESSIVE OFFSET, SPLIT DIMENSIONS
  //cudaMalloc((void **)&m_dH_detJ_dx, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));
  //cudaMalloc((void **)&m_dH_detJ_dy, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));  
  //cudaMalloc((void **)&m_dH_detJ_dz, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));  
  malloc_t (m_H,          double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  
  malloc_t (m_dH_detJ_dx,double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  malloc_t (m_dH_detJ_dy,double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  malloc_t (m_dH_detJ_dz,double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  

  malloc_t(m_detJ,  double, m_elem_count * m_gp_count );    
  malloc_t(m_ematm, double, m_nodxelem * m_nodxelem * m_elem_count); //Elemental mas matrices

  //cudaMalloc((void **)&m_detJ,  m_elem_count * m_gp_count * sizeof (double)); 

  //cudaMalloc((void **)&m_str_rate,  m_elem_count * m_gp_count * 6 * sizeof (double)); 
  //cudaMalloc((void **)&m_rot_rate,  m_elem_count * m_gp_count * 6 * sizeof (double)); 
  //cudaMalloc((void **)&m_sigma,     m_elem_count * m_gp_count * 6 * sizeof (double)); 
  malloc_t(m_str_rate,  double, 6 * m_elem_count * m_gp_count );   
  malloc_t(m_rot_rate,  double, 6 * m_elem_count * m_gp_count );     
  malloc_t(m_sigma,     double, 6 * m_elem_count * m_gp_count );   
  malloc_t(m_tau,       double, 6 * m_elem_count * m_gp_count );   
  
  ///// ELEMENTAL VALUES
  // cudaMalloc((void **)&p,         m_elem_count * m_gp_count * sizeof (double)); 
  // cudaMalloc((void **)&rho,       m_elem_count * m_gp_count * sizeof (double)); 
  // cudaMalloc((void **)&rho_0,     m_elem_count * m_gp_count * sizeof (double)); 
  malloc_t(p,      double, m_elem_count * m_gp_count ); 
  malloc_t(rho,    double, m_elem_count * m_gp_count );   
  malloc_t(rho_0,  double, m_elem_count * m_gp_count ); 
  
  malloc_t(vol,      double, m_elem_count); 
  malloc_t(vol_0,    double, m_elem_count); 
  // cudaMalloc((void **)&vol,       m_elem_count * sizeof (double)); 
  // cudaMalloc((void **)&vol_0,     m_elem_count * sizeof (double)); 

  malloc_t(m_f_elem,    double, m_elem_count * m_dim * m_nodxelem);   
  malloc_t(m_f_elem_hg, double, m_elem_count * m_dim * m_nodxelem);   
  //cudaMalloc((void **)&m_f_elem,  m_elem_count * m_dim * m_nodxelem * sizeof (double)); 

  malloc_t(mat,    Material_*, m_elem_count);   
  //cudaMalloc((void**)&mat,    m_elem_count * sizeof(Material_ *));
  
  #ifdef CUDA_BUILD
	report_gpu_mem_();
  #endif

  malloc_t (x_h,      double,node_count*3);
  //malloc_t (u_h,      double,node_count*3);
  u_h = new double [3*m_node_count];
}


void Domain_d::AssignMaterial (Material_ *material_h) {
//    cudaMalloc((void**)&materials, 1 * sizeof(Material_ )); //
    malloc_t(materials, Material_,1);
    memcpy_t(materials, material_h, 1 * sizeof(Material_));	
}

dev_t void Domain_d::AssignMatAddress(){
  par_loop(i, m_elem_count)
    mat[i] = &materials[0];
  
}

///// MAKE TEMPLATE SET VECTOR
void Domain_d::setDensity(const double &r){
  double *rho_h = new double [m_node_count];
  for (int n=0;n<m_elem_count;n++) rho_h[n] = r;
  memcpy_t(this->rho_0, rho_h, sizeof(double) * m_elem_count);    
  
  delete rho_h;
  
}

dev_t void Domain_d::UpdatePrediction(){
  par_loop (n,m_node_count){
    // vector_t p_a = Ptr_vector_t(prev_a, n);

    // vector_t u_ = dt * (getV(n) + (0.5 - m_beta)* dt *p_a) ;// = dt * (getV(n) + 0.5 - m_beta);

    // vector_t x_ = Ptr_vector_t(x, n);

    // vector_t_Ptr(u_,u_dt,n);
    // vector_t v_ = getV(n) + (1.0 - m_gamma) * dt * p_a; //nod%v = nod%v + (1.0d0-gamma)* dt * prev_a
    // vector_t_Ptr(v_,v,n);

    // PREDICTION PHASE

  }

    for (int i = 0; i < m_node_count; i++) {
        for (int j = 0; j < m_dim; j++) {
            //u_[i][j] = dt * (v_[i][j] + (0.5 - m_beta) * dt * prev_a_[i][j]);
            //NEW; global
            int ig = i*m_dim + j;
            u_dt[ig] = dt * (v[ig] + (0.5 - m_beta) * dt * prev_a[ig]);
        }
    }

    for (int i = 0; i < m_node_count; i++) {
        for (int j = 0; j < m_dim; j++) {
            //v_[i][j] += (1.0 - m_gamma) * dt * prev_a_[i][j];
            
            v[m_dim*i+j] += (1.0 - m_gamma) * dt * prev_a[m_dim*i+j];                
            //printf("v %e",v[m_dim*i+j] );
        }
    }


}

//////////////////////////////////////
////////// THIS NOT INCLUDE DISPLACEMENT BECAUSE OF 
////////// VELOCITY BC ARE PARALLELLIZED BY BCs
////////// INSTEAD OF NODES

dev_t void Domain_d::UpdateCorrectionAccVel(){
  double f = 1.0/(1.0-m_alpha);
  // par_loop (n,m_node_count){

    // vector_t p_a = Ptr_vector_t(prev_a, n);    

    // vector_t a_ = f*(Ptr_vector_t(a, n) - m_alpha * p_a);

    // vector_t_Ptr(a_,a,n);
    // vector_t_Ptr(a_,prev_a,n);
    // vector_t v_ = getV(n) + m_gamma * dt * a_;
    // vector_t_Ptr(v_,v,n);

  // }


    for (int i = 0; i < m_node_count; i++) {
        for (int j = 0; j < m_dim; j++) {
            int ig = i*m_dim + j;

            //printf ("a ig %f\n",  a[ig]);
            a[ig] = f * a[ig]  -m_alpha * prev_a[ig]; //GLOBAL  

            v[ig] += m_gamma * dt * a[ig];
            
        }
    }
        

}


  // !u = u + beta * nod%v * dt
  // u = u + beta * dt * dt * nod%a   
  // nod%u = nod%u + u
  // nod%x = nod%x + u

dev_t void Domain_d   ::UpdateCorrectionPos(){
  double f = 1.0/(1.0-m_alpha);
  // par_loop (n,m_node_count){
    // vector_t uinc_  = Ptr_vector_t(u_dt, n) + m_beta * dt * dt * Ptr_vector_t(prev_a, n); // = dt * (getV(n) + 0.5 - m_beta);
    // vector_t u_     = Ptr_vector_t(u, n) + uinc_;
    // vector_t x_     = Ptr_vector_t(x, n);
    // vector_t_Ptr(u_+x_,x,n);
    // vector_t_Ptr(u_,u,n);       //Copy displacements to device
    

    // vector_t xc_  =Ptr_vector_t(x, n);
  // }

      for (int i = 0; i < m_node_count; i++) {
          for (int j = 0; j < m_dim; j++) {
              // u_[i][j] += m_beta * dt * dt * a_[i][j];
              // x_[i][j] += u_[i][j];
              
              int ig = i*m_dim + j;
              u_dt[ig] += m_beta * dt * dt * a[ig];
              x[ig] += u_dt[ig];
          }
      }

      for (int i = 0; i < m_node_count; i++) {
          for (int j = 0; j < m_dim; j++) {
              //prev_a_[i][j] = a_[i][j];
              
              prev_a[m_dim*i+j] = a[m_dim*i+j];
          }
      }

      for (int i = 0; i < m_node_count; i++) {
          for (int j = 0; j < m_dim; j++) {
              //u_tot_[i][j] += u_[i][j];
              u[m_dim*i+j] += u_dt[m_dim*i+j];
              //printf ("U %.6e \n", u[m_dim*i+j] );
          }
      }
        
}

host_ void Domain_d::ImposeBCAAllDim()//ALL DIM
{    
  for (int d=0;d<m_dim;d++){
    int N = bc_count[d];
    
    #ifdef CUDA_BUILD
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    ImposeBCAKernel<<<blocksPerGrid,threadsPerBlock >>>(this, d);
    cudaDeviceSynchronize();
    #else
      ImposeBCA(d);
    #endif
  }
}

host_ void Domain_d::ImposeBCVAllDim()//ALL DIM
{    
  for (int d=0;d<m_dim;d++){
    int N = bc_count[d];
    
    #ifdef CUDA_BUILD
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    ImposeBCVKernel<<<blocksPerGrid,threadsPerBlock >>>(this, d);
    cudaDeviceSynchronize();
    #else
      ImposeBCV(d);
    #endif
  }
}

host_   void Domain_d::AddBCVelNode(const int &node, const int &dim, const double &val){
  if (dim == 0)       { bcx_nod_h.push_back(node);bcx_val_h.push_back(val);bc_count[0]++;}
  else if (dim == 1)  { bcy_nod_h.push_back(node);bcy_val_h.push_back(val);bc_count[1]++;}
  else if (dim == 2)  { bcz_nod_h.push_back(node);bcz_val_h.push_back(val);bc_count[2]++;}
}

host_   void Domain_d::AllocateBCs() {
  for (int d=0;d<m_dim;d++){cout << "Allocated "<<bc_count[d]<< " Velocity BCs for dof "<<d<<endl; }
  // cudaMalloc((void **)&bcx_nod, bc_count[0] * sizeof (int)); 
  // cudaMalloc((void **)&bcy_nod, bc_count[1] * sizeof (int)); 
  // cudaMalloc((void **)&bcz_nod, bc_count[2] * sizeof (int)); 

  malloc_t(bcx_nod, int, bc_count[0] );     
  malloc_t(bcy_nod, int, bc_count[1] );     
  malloc_t(bcz_nod, int, bc_count[2] );     

  malloc_t(bcx_val, double, bc_count[0] );     
  malloc_t(bcy_val, double, bc_count[1] );     
  malloc_t(bcz_val, double, bc_count[2] );     
  
  // cudaMalloc((void **)&bcx_val, bc_count[0] * sizeof (double)); 
  // cudaMalloc((void **)&bcy_val, bc_count[1] * sizeof (double)); 
  // cudaMalloc((void **)&bcz_val, bc_count[2] * sizeof (double)); 
  
  //ONLY FOR COMPLETENESS IN CASE OF CPU CODE
  double *bcx_val_h_ptr = new double [bc_count[0]];
  int *bcx_nod_h_ptr = new int [bc_count[0]];

  double  *bcy_val_h_ptr = new double [bc_count[1]];
  int     *bcy_nod_h_ptr = new int [bc_count[1]];

  double  *bcz_val_h_ptr = new double [bc_count[2]];
  int     *bcz_nod_h_ptr = new int    [bc_count[2]];  

  for (int i=0;i<bcx_nod_h.size();i++){bcx_nod_h_ptr[i]= bcx_nod_h[i];bcx_val_h_ptr[i]= bcx_val_h[i];}
  for (int i=0;i<bcy_nod_h.size();i++){bcy_nod_h_ptr[i]= bcy_nod_h[i];bcy_val_h_ptr[i]= bcy_val_h[i];} 
  for (int i=0;i<bcz_nod_h.size();i++){bcz_nod_h_ptr[i]= bcz_nod_h[i];bcz_val_h_ptr[i]= bcz_val_h[i];}
  
  memcpy_t(bcx_val, bcx_val_h_ptr, sizeof(double) * bc_count[0]);    
  memcpy_t(bcx_nod, bcx_nod_h_ptr, sizeof(int) * bc_count[0]);   
  memcpy_t(bcy_val, bcy_val_h_ptr, sizeof(double) * bc_count[1]);    
  memcpy_t(bcy_nod, bcy_nod_h_ptr, sizeof(int) * bc_count[1]);   
  memcpy_t(bcz_val, bcz_val_h_ptr, sizeof(double) * bc_count[2]);    
  memcpy_t(bcz_nod, bcz_nod_h_ptr, sizeof(int) * bc_count[2]);   

  delete bcx_val_h_ptr,bcy_val_h_ptr,bcz_val_h_ptr;
  delete bcx_nod_h_ptr,bcy_nod_h_ptr,bcz_nod_h_ptr;
  
}

dev_t void Domain_d::ImposeBCV(const int dim){
  par_loop (n,bc_count[dim]){
    double val;
    //printf("thread %d, Imposing Vel in dim %d, %d Conditions, val %f\n", n, dim, bc_count[dim], bcx_val[n]);
    //printf("BCV dim %d\n", dim);
    if (dim == 0)       {/*printf ("val %f, Nod %d\n",bcx_val[n],bcx_nod[n]); */ v[3*bcx_nod[n]+dim] = bcx_val[n]; }
    else if (dim == 1)  {/*printf ("node %d val %f \n",bcy_nod[n], bcy_val[n]);*/v[3*bcy_nod[n]+dim] = bcy_val[n];}
    else if (dim == 2)  {/*printf ("val %f, Nod %d\n",bcz_val[n],bcz_nod[n]); */ v[3*bcz_nod[n]+dim] = bcz_val[n]; }
  }
  
}

dev_t void Domain_d::ImposeBCA(const int dim){
  par_loop (n,bc_count[dim]){
    double val;
    //printf("thread %d, Imposing Vel in dim %d, %d Conditions\n", n, dim, bc_count[dim]);
    
    if (dim == 0)       {/*//printf ("val %f, Nod %d\n",bcx_val[n],bcx_nod[n]);*/ a[3*bcx_nod[n]+dim] = 0.0; }
    else if (dim == 1)  {/*//printf ("val %f \n",bcy_val[n]);*/                   a[3*bcy_nod[n]+dim] = 0.0;}
    else if (dim == 2)  {/*//printf ("val %f, Nod %d\n",bcz_val[n],bcz_nod[n]); */a[3*bcz_nod[n]+dim] = 0.0; }
  }
  
}

void Domain_d::AddBoxLength(vector_t const & V, vector_t const & L, const double &r,const bool &red_int){
    // integer, intent(in):: tag
  // logical, intent(in) :: redint
  // !real(fp_kind), intent(in), allocatable :: V
  // real(fp_kind), dimension(1:3), intent(in)  :: V ! input
  // real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
  vector_t Xp;
  int p, nnodz;

  int nel[3];

  m_dim = 2;
  if (L.z > 0.0) m_dim = 3;
  
  
  nel[0] = (int)(L.x/(2.0*r));
  nel[1] = (int)(L.y/(2.0*r));
  cout << "Nel x: "<<nel[0]<<", y "<<nel[1]<<endl;
  
  m_gp_count = 1;
  cout << "Mesh dimension is: "<<m_dim<<endl;
  if (m_dim == 2){
    nel[2] = 1;
    m_nodxelem = 4;
    cout << "ne "<<m_nodxelem<<endl;
    if (!red_int) m_gp_count = 4;
  } else {
    nel[2] = (int)(L.z/(2.0*r));
    m_nodxelem = 8;
    if (!red_int) m_gp_count = 8; 
  }
  

  Xp.z = V.z ;
    

    // write (*,*) "Creating Mesh ...", "Elements ", neL.y, ", ",neL.z
    int nc;
  if (m_dim == 2) nc = (nel[0] +1) * (nel[1]+1);
  else            nc = (nel[0] +1) * (nel[1]+1) * (nel[2]+1);
  int ne = nel[0]*nel[1]*nel[2];
  cout << "Mesh created. Element count: "<< nel[0]<<", "<<nel[1]<<", "<<nel[2]<<endl;
  
  //thisAllocateNodes((nel[0] +1) * (nel[1]+1) * (nel[2]+1));
    // print *, "Element count in XYZ: ", nel(:)
    // write (*,*) "Box Node count ", node_count

	this->SetDimension(nc,ne);	 //AFTER CREATING DOMAIN
  cout << "Mesh generated. Node count: " << nc<<". Element count: "<<ne<<endl;
  cout << "Dimension is: "<<m_dim<<endl;
  //SPH::Domain	dom;
	//vector_t *x =  (vector_t *)malloc(dom.Particles.size());
	vector_t *x_H =  new vector_t [m_node_count];


	//int size = dom.Particles.size() * sizeof(vector_t);
	cout << "Copying to device..."<<endl;
    
  cout << "Box Particle Count is " << m_node_count <<endl;
  p = 0;
  int kmax = nel[2] +1;
  if (m_dim == 2) kmax = 1;
  for (int k = 0; k < kmax;k++) {
    Xp.y = V.y;
    for (int j = 0; j < (nel[1] +1);j++){
      Xp.x = V.x;
      for (int i = 0; i < (nel[0] +1);i++){
        //m_node.push_back(new Node(Xp));
        x_H[p] = Xp;
        //nod%x(p,:) = Xp(:);
        cout << "node " << p <<"X: "<<Xp.x<<"Y: "<<Xp.y<<"Z: "<<Xp.z<<endl;
        p++;
        Xp.x = Xp.x + 2.0 * r;
      }
      Xp.y = Xp.y + 2.0 * r;
    }// 
    Xp.z = Xp.z + 2 * r;

  //cout <<"m_node size"<<m_node.size()<<endl;
  } 
		memcpy_t(this->x,   x_H, sizeof(vector_t) * m_node_count);    

    // !! ALLOCATE ELEMENTS
    // !! DIMENSION = 2
    int gp = 1;
    if (m_dim == 2) {
      // if (redint .eqv. .False.) then
        // gp = 4
      // end if 
      //call AllocateElements(neL.y * neL.z,gp) !!!!REDUCED INTEGRATION
    } else {
      // if (redint .eqv. .False.) then
        // gp = 8
      // end if 
      // call AllocateElements(neL.y * neL.z*nel(3),gp) 
    }

		elnod_h       = new int [m_elem_count * m_nodxelem]; //Flattened

     int *nodel_count_h  = new int [m_node_count];
     int *nodel_offset_h = new int [m_node_count];
    for (int n=0;n<m_node_count;n++){
      nodel_count_h[n] = 0;
    }
    
    cout << "Allocating element nodes.."<<endl;
		int ex, ey, ez;
		std::vector <int> n;
    if (m_dim == 2) {
			n.resize(4);
      int ei = 0;
      for (int ey = 0; ey < nel[1];ey++){
        for (int ex = 0; ex < nel[0];ex++){
        int iv[4];
        int nb1 = (nel[0]+1)* ey    + ex;    
        int nb2 = (nel[0]+1)*(ey+1) + ex;
        elnod_h[ei  ] = nb1;                        nodel_count_h[nb1  ] ++;             
        elnod_h[ei+1] = nb1 + 1;                    nodel_count_h[nb1+1] ++;             
        elnod_h[ei+2] = nb2 + 1;                    nodel_count_h[nb2+1] ++;      
        elnod_h[ei+3] = (nel[0]+1)*(ey+1) + ex;     nodel_count_h[nb2 ] ++;      
			
				 for (int i=0;i<m_nodxelem;i++)cout << elnod_h[ei+i]<<", ";
					cout << "Nel x : "<<nel[0]<<endl;
					cout << "nodes "<<endl;
					ei += m_nodxelem;
					 }
      } 
    } else { //dim: 3
      int ei = 0; //ELEMENT INTERNAL NODE (GLOBAL INDEX)
      int nnodz = (nel[0]+1)*(nel[1]+1);
      for (int ez = 0; ez < nel[2];ez++)
      for (int ey = 0; ey < nel[1];ey++){
        for (int ex = 0; ex < nel[0];ex++){
          
          int iv[8];
          int nb1 = nnodz*ez + (nel[0]+1)*ey + ex;
          int nb2 = nnodz*ez + (nel[0]+1)*(ey+1) + ex;
          
          elnod_h[ei  ] = nb1;                      nodel_count_h[nb1  ] ++;          
          elnod_h[ei+1] = nb1+1;                    nodel_count_h[nb1+1] ++;
          elnod_h[ei+2] = nb2+1;                    nodel_count_h[nb2+1] ++;
          elnod_h[ei+3] = nb2;                      nodel_count_h[nb2  ] ++;
          
          elnod_h[ei+4] = nb1 + nnodz;              nodel_count_h[nb1 + nnodz    ]++;   
          elnod_h[ei+5] = nb1 + nnodz + 1;          nodel_count_h[nb1 + nnodz + 1]++;  
          elnod_h[ei+6] = nb2 + nnodz + 1;          nodel_count_h[nb2 + nnodz + 1]++;  
          elnod_h[ei+7] = nb2 + nnodz;              nodel_count_h[nb2 + nnodz    ]++;  
          
          // for (int i=0;i<8;i++)
            // cout << elnod_h[ei + i]<<", ";
          // cout <<endl;

            cout << "Nel x : "<<nel[0]<<", Element: " << ei/m_nodxelem<<endl;
           cout << "nodes "<<endl;
           cout << "nodxelem " <<m_nodxelem<<endl;
           cout << "nb1, nb2 "<< nb1 << ", "<<nb2<<" nnodz"<< nnodz*(ez+1)<<"ez "<<ez<<endl;
           for (int i=0;i<m_nodxelem;i++)cout << elnod_h[ei+i]<<", ";
           ei += m_nodxelem;

					 }
      } 

		}//if dim 

    //cudaMalloc((void **)&m_elnod, m_elem_count * m_nodxelem * sizeof (int));	
    malloc_t(m_elnod, unsigned int, m_elem_count * m_nodxelem);
		memcpy_t(this->m_elnod, elnod_h, sizeof(int) * m_elem_count * m_nodxelem); 
    
    // printf("ELNOD \n");
    // for (int e=0;e<m_elem_count;e++){
      // for (int n=0;n<m_nodxelem;n++)
        // printf("%d ",m_elnod[e*m_nodxelem+n]);
      // printf("\n");
    // }
    //cudaMalloc(&m_jacob,m_elem_count * sizeof(Matrix ));
    malloc_t(m_jacob, Matrix, m_elem_count );
    
    //////////////////// ELEMENT SHARED BY NODES (FOR PARALLEL NODAL MODE ASSEMBLY) ///////////////////////////////
    int nodel_tot = 0;
    for (int n=0;n<m_node_count;n++){
      nodel_offset_h[n] = nodel_tot;
      nodel_tot        += nodel_count_h[n];
      cout << "NodEL tot " << nodel_tot<<endl;
      cout << "Node "<< n << " Shared elements: "<<nodel_count_h[n]<<endl;

    }
    cout << "Size of Nodal shared Elements vector "<< nodel_tot<<endl;
		int *nodel_h       = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT
    int *nodel_loc_h   = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT    
    
    for (int n=0;n<m_node_count;n++)  nodel_count_h[n] = 0;    
    for (int e=0;e<m_elem_count;e++){
      int offset = m_nodxelem * e;
      for (int ne=0;ne<m_nodxelem;ne++){
        int n = elnod_h[offset+ne];
        
        nodel_h     [nodel_offset_h[n] + nodel_count_h[n]] = e;
        nodel_loc_h [nodel_offset_h[n] + nodel_count_h[n]] = ne;
        
        nodel_count_h[n]++;
      }//nod x elem 
    }

    // cudaMalloc((void **)&m_nodel,     nodel_tot * sizeof (int));
    // cudaMalloc((void **)&m_nodel_loc, nodel_tot * sizeof (int));
    
    malloc_t (m_nodel,        int,nodel_tot);
    malloc_t (m_nodel_loc,    int,nodel_tot);
    
    malloc_t (m_nodel_offset, int, m_node_count);
    malloc_t (m_nodel_count,  int, m_node_count);
    
    //THIS IS ONLY FOR COMPLETENESS IN CASE OF CPU, SHOULD BE BETTER TO WRITE ON FINALL ARRAY
		memcpy_t(this->m_nodel,         nodel_h,        sizeof(int) * nodel_tot); 
		memcpy_t(this->m_nodel_loc,     nodel_loc_h,    sizeof(int) * nodel_tot);
		 
		memcpy_t(this->m_nodel_offset,  nodel_offset_h, sizeof(int) * m_node_count);  //OFFSET FOR PREVIOUS ARRAYS
		memcpy_t(this->m_nodel_count,    nodel_count_h, sizeof(int) * m_node_count);  //OFFSET FOR PREVIOUS ARRAYS
		    
    // ///// TESTING
    for (int n=0;n<m_node_count;n++){
      cout << "M node offset:"<<nodel_offset_h[n];
      cout << "Node  "<< n << " Elements"<<endl;
      for (int ne=0;ne<nodel_count_h[n];ne++) cout << nodel_h[nodel_offset_h[n]+ne]<<", ";
      cout << endl;
      cout << "Node  "<< n << " Elements Internal Node"<<endl;
      for (int ne=0;ne<nodel_count_h[n];ne++) cout << nodel_loc_h[nodel_offset_h[n]+ne]<<", ";
      cout << endl;
    }
    
    cout << "Mesh generation done. "<<endl;

		
		delete [] /*elnod_h, */nodel_count_h, nodel_h, nodel_loc_h,nodel_offset_h;
}
///////////////////////////////
//// SINGLE ELNOD AND NODEL CONNECTIVITY
/// ASSUMES 
// void Domain_d::SetConnectivity(int elnod_h[], int elsize, int nodxelem){
  // m_nodxelem = nodxelem;
    // //////////////////// ELEMENT SHARED BY NODES (FOR PARALLEL NODAL MODE ASSEMBLY) ///////////////////////////////
  // int nodel_tot = 0;
  // for (int n=0;n<m_node_count;n++){
    // nodel_offset_h[n] = nodel_tot;
    // nodel_tot        += nodel_count_h[n];
    // cout << "NodEL tot " << nodel_tot<<endl;
    // cout << "Node "<< n << " Shared elements: "<<nodel_count_h[n]<<endl;

  // }
  // cout << "Size of Nodal shared Elements vector "<< nodel_tot<<endl;
  // int *nodel_h       = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT
  // int *nodel_loc_h   = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT    
  
  // for (int n=0;n<m_node_count;n++)  nodel_count_h[n] = 0;    
  // for (int e=0;e<m_elem_count;e++){
    // int offset = m_nodxelem * e;
    // for (int ne=0;ne<m_nodxelem;ne++){
      // int n = elnod_h[offset+ne];
      
      // nodel_h     [nodel_offset_h[n] + nodel_count_h[n]] = e;
      // nodel_loc_h [nodel_offset_h[n] + nodel_count_h[n]] = ne;
      
      // nodel_count_h[n]++;
    // }//nod x elem 
  // }
// }


void Domain_d::CreateFromLSDyna(const lsdynaReader &reader){

  this->SetDimension(reader.m_node.size(),reader.m_elem_count);	 //AFTER CREATING DOMAIN
  cout << "Node Size: "<<m_node_count<<endl;  
  cout << "Element Size: "<<m_elem_count<<endl;  
  
}

dev_t void Domain_d::calcElemJAndDerivatives () {
	////printf("calculating\n");
	////printf ("threadIdx.x %d, blockDim.x%d, blockIdx.x %d\n",threadIdx.x ,blockDim.x , blockIdx.x);
  //printf("allocating m_dim %d, n_nodxelem %d\n", m_dim, m_nodxelem);
  
	Matrix *jacob = new Matrix(m_dim, m_dim);    
	Matrix *inv_j = new Matrix(m_dim, m_dim);    
  //Matrix *dHrs = new Matrix(m_dim, m_nodxelem);   /////////////////////////////// IF CREATION IS DYNAMIC ! (TEST IF )
  
  Matrix *dHxy_detJ_loc = new Matrix(m_dim, m_nodxelem);
   
  //printf("done\n");
   
	////printf ("e %d, elem_count %d\n",m_elem_count);
  par_loop (e, m_elem_count) {
  int offset = m_gp_count * e;
  Matrix *x2 = new Matrix(m_nodxelem, m_dim);
  // integer :: e
  // ! !rg=gauss[ig]
  // ! !sg=gauss[jg]
  // real(fp_kind), dimension(dim,m_nodxelem) :: dHrs !!! USED ONLY FOR SEVERAL GAUSS POINTS
	//printf ("m_dim %d, nod x elem %d", m_dim, m_nodxelem);

  
  //cudaMalloc((void**)&dHrs_p, sizeof(Matrix));
	////printf("test %lf",dHrs.m_data[0]);
	//double dHrs_fl[m_dim* m_nodxelem];
	//dHrs->Print();

  //printf("x2 dimensions %d X %d\n", m_nodxelem, m_dim);
   
   // //printf("Jacob\n");jacob->Print();
   double gpc[8][3];

	//printf ("Matrices created\n");

      // do i=1,nodxelem
          // !print *, "elnod " , elem%elnod(e,i)
          // x2(i,:)=nod%x(elem%elnod(e,i),:)
      // end do
  int nind = e * m_nodxelem;
  for (int i=0;i<m_nodxelem;i++){
      vector_t x_ = Ptr_vector_t(x,m_elnod[nind+i]);
      //printf("elnod %d\n",m_elnod[nind+i]);
      //printf("x: %.6e %.6e %.6e\n",x_.x, x_.y,x_.z);
      x2->Set(i,0,x_.x); x2->Set(i,1,x_.y); 
      if (m_dim == 3)
        x2->Set(i,2,x_.z);
      
      ////printf ("elnod %d, %lf %lf %lf \n",m_elnod[nind+i],x[m_elnod[nind+i]].x,x[m_elnod[nind+i]].y,x[m_elnod[nind+i]].z);
  } 
  //printf("x2\n");x2->Print();
  //printf("m_gp_count %d\n",m_gp_count);
    //printf("Calculating jacobian\n");
    if (m_gp_count == 1 ) {      
      //invJ = adj(elem%jacob(e,gp,:,:)) !!! IN FACT IS invJ x detJ
			if (m_dim == 2) {
      // if (dim .eq. 2) then 
        // !dHdrs [-1,1,1,-1;  -1.-1,1,1] x X2
        // !! J = [
        // !! dx/dr dy/dr
        // !! dx/ds dy/dx ]
        // !!! THIS IS TO AVOID MATMUL
        // ! print *, "nodes X ", x2(:,1)
        // ! print *, "nodes Y ", x2(:,2)
        for (int d=0;d<2;d++){
          // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)
          // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)
          // elem%jacob(e,gp,:,:) = 0.25*elem%jacob(e,gp,:,:)
          jacob->Set(0,d,0.25*(-x2->getVal(0,d)+x2->getVal(1,d)+x2->getVal(2,d)-x2->getVal(3,d))); 
          jacob->Set(1,d,0.25*(-x2->getVal(0,d)-x2->getVal(1,d)+x2->getVal(2,d)+x2->getVal(3,d)));
        }
        
        AdjMat(*jacob, inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
        //printf(" J ptr\n");
        //jacob->Print();
        //printf("ADJ J ptr\n");
        //inv_j->Print();          //printf("jacob\n");jacob->Print();
        for (int d=0;d<2;d++){        
          dHxy_detJ_loc->Set(d,0,0.25*(-inv_j->getVal(d,0)-inv_j->getVal(d,1)));     
          dHxy_detJ_loc->Set(d,1,0.25*(inv_j->getVal(d,0)-inv_j->getVal(d,1)));     
          dHxy_detJ_loc->Set(d,2,0.25*( inv_j->getVal(d,0)+inv_j->getVal(d,1)));     
          dHxy_detJ_loc->Set(d,3,0.25*(-inv_j->getVal(d,0)+inv_j->getVal(d,1)));     
          // #ifdef CUDA_BUILD
          // jacob->Set(0,d,0.25*(-x2->getVal(0,d) + x2->getVal(1,d) + x2->getVal(2,d) - x2->getVal(3,d)));  
          // jacob->Set(1,d,0.25*(-x2->getVal(0,d) - x2->getVal(1,d) + x2->getVal(2,d) + x2->getVal(3,d)));  
          // #else 
          // jacob->operator()(0,d) = -x2->getVal(0,d) + x2->getVal(1,d) + x2->getVal(2,d) - x2->getVal(3,d);
          // #endif

          // elem%dHxy_detJ(e,gp,:,1) = -invJ(:,1)-invJ(:,2) !For each 3 rows of inv J and dHdxy
          // elem%dHxy_detJ(e,gp,:,2) =  invJ(:,1)-invJ(:,2)
          // elem%dHxy_detJ(e,gp,:,3) =  invJ(:,1)+invJ(:,2)
          // elem%dHxy_detJ(e,gp,:,4) = -invJ(:,1)+invJ(:,2)     
          
          // elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.25d0
        }
        //dHxy_detJ_loc->Mul(0.25);
        
			} else { //!!!DIM 3

          for (int d=0;d<m_dim;d++){
            jacob->Set(0,d,0.125*(-x2->getVal(0,d)+x2->getVal(1,d)+x2->getVal(2,d)-x2->getVal(3,d)-x2->getVal(4,d)+x2->getVal(5,d)+x2->getVal(6,d)-x2->getVal(7,d)));  
            jacob->Set(1,d,0.125*(-x2->getVal(0,d)-x2->getVal(1,d)+x2->getVal(2,d)+x2->getVal(3,d)-x2->getVal(4,d)-x2->getVal(5,d)+x2->getVal(6,d)+x2->getVal(7,d)));  
            jacob->Set(2,d,0.125*(-x2->getVal(0,d)-x2->getVal(1,d)-x2->getVal(2,d)-x2->getVal(3,d)+x2->getVal(4,d)+x2->getVal(5,d)+x2->getVal(6,d)+x2->getVal(7,d))); 
            //jacob->Set(0,d,-x2->getVal(0,d) + x2->getVal(1,d) + x2->getVal(2,d) - x2->getVal(3,d));  

          }
          //printf("jacob\n");jacob->Print();
          
          // //printf ("Setting dhxy\n");
          //*inv_j = jacob->Inv();
          //Matrix inv = jacob->Inv();
          //InvMat(*jacob, inv_j);
          AdjMat(*jacob, inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
          //printf("ADJ J ptr\n");
          //inv_j->Print();          //printf("jacob\n");jacob->Print();
                  
          // jacob->Print();
          ////printf("INV J2 not ptr\n");
          //inv.Print();
          
          //inv.Print();
          for (int d=0;d<m_dim;d++){            
            dHxy_detJ_loc->Set(d,0,0.125*(-inv_j->getVal(d,0)-inv_j->getVal(d,1)-inv_j->getVal(d,2)));         
            dHxy_detJ_loc->Set(d,1,0.125*( inv_j->getVal(d,0)-inv_j->getVal(d,1)-inv_j->getVal(d,2)));  
            dHxy_detJ_loc->Set(d,2,0.125*( inv_j->getVal(d,0)+inv_j->getVal(d,1)-inv_j->getVal(d,2)));  
            dHxy_detJ_loc->Set(d,3,0.125*(-inv_j->getVal(d,0)+inv_j->getVal(d,1)-inv_j->getVal(d,2)));             
            dHxy_detJ_loc->Set(d,4,0.125*(-inv_j->getVal(d,0)-inv_j->getVal(d,1)+inv_j->getVal(d,2))); 
            dHxy_detJ_loc->Set(d,5,0.125*( inv_j->getVal(d,0)-inv_j->getVal(d,1)+inv_j->getVal(d,2)));
            dHxy_detJ_loc->Set(d,6,0.125*( inv_j->getVal(d,0)+inv_j->getVal(d,1)+inv_j->getVal(d,2)));
            dHxy_detJ_loc->Set(d,7,0.125*(-inv_j->getVal(d,0)+inv_j->getVal(d,1)+inv_j->getVal(d,2)));
          }
          //dHxy_detJ_loc->Mul(0.125); /////->DO NOT USE THIS!! --- ERRORS ---

          // // elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.125d0    
          
      } // end if  !!!!DIM
      
      m_detJ[offset] = jacob->calcDet();
      //printf("det J %f allocated, offset %d\n",m_detJ[offset],offset);
      // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
    } else { //!!!!! GP > 1
			
      // double r = 1.0/sqrt(3.0);
			// gpc[0][0] = -r; gpc[0][1] = -r;gpc[0][2] = -r;
			// gpc[1][0] =  r; gpc[1][1] = -r;gpc[1][2] = -r;
			// gpc[2][0] = -r; gpc[2][1] =  r;gpc[2][2] = -r;
			// gpc[3][0] =  r; gpc[3][1] =  r;gpc[3][2] = -r;
			// gpc[4][0] = -r; gpc[4][1] = -r;gpc[4][2] =  r;
			// gpc[5][0] =  r; gpc[5][1] = -r;gpc[5][2] =  r;
			// gpc[6][0] = -r; gpc[6][1] =  r;gpc[6][2] =  r;
			// gpc[7][0] =  r; gpc[7][1] =  r;gpc[7][2] =  r;
			
			// //,:)=[-r,-r,-r];   gpc(2,:)=[ r,-r,-r];      gpc(3,:)=[-r, r,-r];      gpc(4,:)=[ r, r,-r]; !These are the 4 points for 2D full elem
      // // gpc(1,:)=[-r,-r,-r];   gpc(2,:)=[ r,-r,-r];      gpc(3,:)=[-r, r,-r];      gpc(4,:)=[ r, r,-r]; !These are the 4 points for 2D full elem
      // // gpc(5,:)=[-r,-r, r];   gpc(6,:)=[ r,-r, r];      gpc(7,:)=[-r, r, r];      gpc(8,:)=[ r, r, r];
      // //h1 = (1-r)(1-s)(1-t) //h2 = (1+r)(1-s)(1-t) //h3 = (1-r)(1+s)
      // //h3 = (1+r)(1+s)(1-t) //h4 = (1-r)(1+s)(1-t)
            // // elem%math(e,gp, 1,:) = 0.125*[(1-gpc(gp,1))*(1-gpc(gp,2))*(1-gpc(gp,3)),(1+gpc(gp,1))*(1-gpc(gp,2))*(1-gpc(gp,3)), &
                                // // (1+gpc(gp,1))*(1+gpc(gp,2))*(1-gpc(gp,3)),(1-gpc(gp,1))*(1+gpc(gp,2))*(1-gpc(gp,3)), &
                                // // (1-gpc(gp,1))*(1-gpc(gp,2))*(1+gpc(gp,3)),(1+gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3)), &
                                // // (1+gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3)),(1-gpc(gp,1))*(1+gpc(gp,2))*(1+gpc(gp,3))]
      // if (m_dim == 3) {
        // for (int gp=0;gp<m_gp_count;gp++){
          
          // dHrs->Set(0,0,-1.0*(1-gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs->Set(1,0,-1.0*(1+gpc[gp][0])*(1.0-gpc[gp][2])); dHrs->Set(2,0,-1.0*(1+gpc[gp][0])*(1.0-gpc[gp][1])); //dh1/d(r,s,t)
          // dHrs->Set(0,1,     (1-gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs->Set(1,1,-1.0*(1+gpc[gp][0])*(1.0-gpc[gp][2])); dHrs->Set(2,1,-1.0*(1-gpc[gp][0])*(1.0-gpc[gp][1])); //dh2/d(r,s,t)
					
					// dHrs->Set(0,2,     (1+gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs->Set(1,2,     (1+gpc[gp][0])*(1.0-gpc[gp][2])); dHrs->Set(2,2,-1.0*(1+gpc[gp][0])*(1.0+gpc[gp][1]));
					// dHrs->Set(0,3,-1.0*(1+gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs->Set(1,3,     (1-gpc[gp][0])*(1.0-gpc[gp][2])); dHrs->Set(2,3,-1.0*(1+gpc[gp][0])*(1.0+gpc[gp][1]));
          
          // dHrs->Set(0,4,-1.0*(1-gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs->Set(1,4,-1.0*(1-gpc[gp][0])*(1.0+gpc[gp][2])); dHrs->Set(2,4,     (1-gpc[gp][0])*(1.0-gpc[gp][1]));
          // dHrs->Set(0,5,     (1-gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs->Set(1,5,-1.0*(1+gpc[gp][0])*(1.0+gpc[gp][2])); dHrs->Set(2,5,     (1+gpc[gp][0])*(1.0-gpc[gp][1]));
          
          // dHrs->Set(0,6,     (1+gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs->Set(1,6,     (1+gpc[gp][0])*(1.0+gpc[gp][2])); dHrs->Set(2,6,     (1+gpc[gp][0])*(1.0+gpc[gp][1]));
          // dHrs->Set(0,7,-1.0*(1+gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs->Set(1,7,     (1-gpc[gp][0])*(1.0+gpc[gp][2])); dHrs->Set(2,7,     (1-gpc[gp][0])*(1.0+gpc[gp][1]));
					          

					// //*jacob = 0.125 * MatMul(*dHrs,*x2);
          // MatMul(*dHrs,*x2,jacob);
          // //printf("x2\n");
          // //m_jacob[e].Print();

          // x2->Print();
          // jacob->Mul(0.125);
          //printf("jacob\n");jacob->Print();
          //jacob->Print();
          // ////printf("Jacobian: \n");jacob->Print();
          // //printf("dHrs\n"); dHrs->Print();
           
          // InvMat(*jacob, inv_j);
          // //printf("inv j\n");inv_j->Print();
          // MatMul(*inv_j,*dHrs,dHxy_detJ_loc);
          
          // //printf("Derivative matrix\n");
          // dHxy_detJ_loc->Print();
         
          
          // m_detJ[offset + m_gp_count * gp] = jacob->calcDet();
          // //printf("det J %f\n",m_detJ[offset + m_gp_count * gp]);
          // //TRY WITHOUT ALLOCATING
          
          // // invJ = adj(elem%jacob(e,gp,:,:))!!!!/elem%detJ(e,gp) !!!! ALREADY CALCULATED    
          // // !print *, "detJ", elem%detJ(e,gp)
          // // !print *, "invJ", invJ
          // // elem%dHxy_detJ(e,gp,:,:) = 0.125d0 * matmul(invJ,elem%dHrs(e,gp,:,:))
          
        // }// gp
      // } else { //!dim =2
        // // do gp = 1,4
          // // dHrs(1,:)=[-1.0*(1-gpc(gp,2)),     (1-gpc(gp,2))&
                    // // ,     (1+gpc(gp,2)),-1.0*(1+gpc(gp,2))]
          // // dHrs(2,:)=[-1.0*(1-gpc(gp,1)),-1.0*(1+gpc(gp,1))&
                         // // ,(1+gpc(gp,1)),     (1-gpc(gp,1))]  
					// for (int gp=0;gp<m_gp_count;gp++){										
						// dHrs->Set(0,0,-1.0*(1-gpc[gp][1])); dHrs->Set(0,1,     (1-gpc[gp][1])); dHrs->Set(0,2, 1+gpc[gp][1]);   dHrs->Set(0,3,-1.0*(1+gpc[gp][1]));
						// dHrs->Set(1,0,-1.0*(1-gpc[gp][0])); dHrs->Set(1,1,-1.0*(1+gpc[gp][0])); dHrs->Set(1,2,(1+gpc[gp][0]));  dHrs->Set(1,3,     (1-gpc[gp][0]));
					// }
					// //*jacob = 0.125 * MatMul(*dHrs,*x2);
          // MatMul(*dHrs,*x2,jacob);
          // //printf("jacob\n");
          // //m_jacob[e].Print();

          // x2->Print();
          // jacob->Print();
          // jacob->Mul(0.125);
          
          // //m_detJ[offset + m_gp_count * gp] = det(*jacob);
          
          // // elem%dHrs(e,gp,:,:) =  dHrs(:,:)         
          // // !dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]         
          // // !dHrs(3,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))] 
          // // !print *, "dhrs", dHrs 
          // // !print *, "x2", x2 
          // // elem%jacob(e,gp,:,:) = 0.25*matmul(dHrs,x2)
					// //*jacob = 0.25 * MatMul(*dHrs,*x2);
					
					
					// //jacob->Print();
        
        
        
      // }//dim 2 (gp>1)
    }// end if !!gp ==1

    ///// ALLOCATION
    for (int gp=0;gp<m_gp_count;gp++){
      //Domain_d::setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &v)
      //setDerivative(e,gp,dHxy_detJ_loc
      for (int j=0;j<m_nodxelem;j++){
        int offset = e*(m_nodxelem * m_gp_count) + gp * m_nodxelem;
        ////printf ("Offset %d \n", offset);
        
          //m_dH_detJ_dx[offset + j                 ] = dHxy_detJ_loc->operator()(0,j);
          // m_dH_detJ_dx[offset + j] = dHxy_detJ_loc->getVal(0,j);
          // m_dH_detJ_dy[offset + j] = dHxy_detJ_loc->getVal(1,j); 
          // m_dH_detJ_dz[offset + j] = dHxy_detJ_loc->getVal(2,j);      
          setDerivative(e,gp,0,j,dHxy_detJ_loc->getVal(0,j));
          setDerivative(e,gp,1,j,dHxy_detJ_loc->getVal(1,j));
          if (m_dim ==3)
            setDerivative(e,gp,2,j,dHxy_detJ_loc->getVal(2,j));
          //printf("set der: z n %d %f\n",j, dHxy_detJ_loc->getVal(2,j));
          
      }
    }
          
    //printf("jacob\n");
    //jacob->Print();
    //printf("dHdx x detJ\n");
    //dHxy_detJ_loc->Print();
		//printf("END.\n");
    delete x2; //DEFINED ON EACH BLOCK!
  } // e < elem_colunt
  
      delete inv_j, jacob,dHxy_detJ_loc;
}

// __device__ double & Domain_d::getDerivative(const int &e, const int &gp, const int &i, const int &j){
  // //int offset = m_nodxelem * m_gp_count;
  // //if (e < m_elem_count) {
      // return m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i];
  // //}
  // //return ret;
// }

int Domain_d::WriteToCSV(char *FileKey){
  std::ostringstream oss;
	//Writing in a Log file
	//String fn(FileKey);
	std::string fn(FileKey);
	
	oss << "X, Y, Z"<<endl;;
  
  double *xh;
  
  xh = x; //IF !CUDA
  
	
	//#pragma omp parallel for schedule(static) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	// #else
	//printf("Dimension: %d",Dimension);
	for (int i=0; i<m_node_count; i++)//Like in Domain::Move
	{
		//for (int j=0;j<3;j++)
      printf ("%.6e ,%.6e ,%.6e ",xh[i*3],xh[i*3+1],xh[i*3+2]);
			oss << xh[i*3]<< ", "<<xh[i*3+1]<<xh[i*3+2]<<endl;
		
		//Particles[i]->CalculateEquivalentStress();		//If XML output is active this is calculated twice
		//oss << Particles[i]->Sigma_eq<< ", "<< Particles[i]->pl_strain <<endl;
	}

	// fn = FileKey;
	// fn.append(".csv");	
	std::ofstream of(fn.c_str(), std::ios::out);
	of << oss.str();
	of.close();
}

dev_t void Domain_d::printVec(double *v){
  for (int n=0;n<m_node_count;n++){
    for (int d=0;d<m_dim;d++)
      printf("%.6e ",v[m_dim*n + d]);
    printf("\n");
  }
}
//ONLY FOR 3D
dev_t void Domain_d::printSymmTens(double *v){
  for (int e=0;e<m_elem_count;e++){
  tensor3 ss;
  for (int gp=0;gp<m_gp_count;gp++){
      int offset_s = e * m_gp_count + gp;   //SCALAR
      int offset_t = offset_s * 6 ; //SYM TENSOR
      ss = FromFlatSym(v,          offset_t );
      print(ss);
  // printf("TENST STRESSES\n");
  // print(ss);
  // for (int n=0;n<m_elem_count;n++){
    // printf ("DIMENSION %d\n", n);
      // printf("%.6e %.6e %.6e\n",v[6*n  ],v[6*n+3],v[6*n+5]);
      // printf("%.6e %.6e %.6e\n",v[6*n+3],v[6*n+1],v[6*n+4]);
      // printf("%.6e %.6e %.6e\n",v[6*n+5],v[6*n+4],v[6*n+2]);      
    // printf("\n");
    }
  }
}


__global__ void calcElemJAndDerivKernel(Domain_d *dom_d){
		
		dom_d->calcElemJAndDerivatives();
}

__global__ void ImposeBCVKernel(Domain_d *dom_d, const int d){
    dom_d->ImposeBCV(d);
}

__global__ void ImposeBCAKernel(Domain_d *dom_d, const int d){
    dom_d->ImposeBCA(d);
}

__global__ void UpdatePredictionKernel(Domain_d *dom_d){
  dom_d->UpdatePrediction();
}

__global__ void UpdateCorrectionAccVelKernel(Domain_d *dom_d){
  dom_d->UpdateCorrectionAccVel();
}

__global__ void UpdateCorrectionPosKernel(Domain_d *dom_d){
  dom_d->UpdateCorrectionPos();
}

__global__ void AssignMatAddressKernel(Domain_d *dom){
  dom->AssignMatAddress();
}

__global__ void printVecKernel(Domain_d *dom_d, double *v){
  dom_d->printVec(v);
}

__global__ void printSymmTensKernel(Domain_d *dom_d, double *v){
  dom_d->printSymmTens(v);
}


};
	
