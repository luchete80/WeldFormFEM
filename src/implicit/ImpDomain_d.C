#include "ImpDomain_d.h"
#include <iostream>
#include <vector>

#include "Matrix.h"
#include <sstream>
#include <fstream> 
#include <iostream>


#include "Tensor3.C"
#include "lsdynaReader.h"


using namespace std;
using namespace LS_Dyna;

namespace MetFEM {
  
#define ELNOD  4   //ORIGINALLY 8
#define FACENOD 3  //ORIGINALLY 4
#define ELFAC  4   //ORIGINALLY 6
//OLD FOT HEXA, CHANGE IT



void ImpDomain_d::SetDimension(const int &node_count, const int &elem_count){
  
  m_node_count = node_count;
  m_elem_count = elem_count;
  
  // NODAL VARIABLES
  
  malloc_t (x,      double,node_count*m_dim);
  malloc_t (v,      double,node_count*m_dim);
  malloc_t (a,      double,node_count*m_dim);
  malloc_t (u,      double,node_count*m_dim);
  malloc_t (u_dt,   double,node_count*m_dim);
  
  malloc_t (prev_a, double,node_count*m_dim);  
	//cudaMalloc((void **)&m_f, node_count * sizeof (double) * 3);
  malloc_t (m_fi,double,node_count*m_dim); //Internal forces
  malloc_t (m_fe,double,node_count*m_dim);
  
  malloc_t (m_mdiag, double,node_count);
  malloc_t (m_mglob, double,node_count*node_count); //TODO: MAKE SPARSE. DEALLOCATED AFER DIAG CALCULATION

  //////if thermal
  
  malloc_t (T,      double,node_count);
  malloc_t(m_dTedt,    double, m_elem_count * m_dim * m_nodxelem);

  /// MATRICES ///
  /// dHxy_detJ: DIM X NODXELEM
  /// TO AVOID EXCESSIVE OFFSET, SPLIT DIMENSIONS
  //cudaMalloc((void **)&m_dH_detJ_dx, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));
  //cudaMalloc((void **)&m_dH_detJ_dy, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));  
  //cudaMalloc((void **)&m_dH_detJ_dz, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));  
  malloc_t (m_H,          double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  
  cout <<"setting deriv allocation size of "<<m_dim * m_nodxelem * m_elem_count * m_gp_count<<endl;
  cout <<"mdim: "<<m_dim << " "<<", Nod x Elem "<<m_nodxelem<<" "<< m_elem_count <<" "<<m_gp_count<<endl;
  malloc_t (m_dH_detJ_dx,double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  malloc_t (m_dH_detJ_dy,double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  malloc_t (m_dH_detJ_dz,double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  

  malloc_t(m_detJ,  double, m_elem_count * m_gp_count );    
  
  malloc_t(m_nodxelem_e,  int, m_elem_count);
  
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
  
  // USEFUL FOR TETRA Average Nodal Pressure 
  malloc_t(m_voln,      double, node_count); 
  malloc_t(m_vol_0n,    double, node_count); 

  malloc_t(m_f_elem,    double, m_elem_count * m_dim * m_nodxelem);   
  malloc_t(m_f_elem_hg, double, m_elem_count * m_dim * m_nodxelem);   
  //cudaMalloc((void **)&m_f_elem,  m_elem_count * m_dim * m_nodxelem * sizeof (double)); 

  malloc_t(mat,    Material_*, m_elem_count);   
  //cudaMalloc((void**)&mat,    m_elem_count * sizeof(Material_ *));
  
  // AXISYMM
  malloc_t (m_radius, double, m_elem_count * m_gp_count);
  
  
  malloc_t(m_jacob, Matrix, m_elem_count );
    
  #ifdef CUDA_BUILD
	report_gpu_mem_();
  #endif

  x_h = new double [m_dim*m_node_count];
  u_h = new double [m_dim*m_node_count];
  
  ////// CONTACT//////////////////////////////////////////////
  malloc_t(ext_nodes, bool, m_node_count);
  malloc_t(contforce, double,node_count*m_dim);
  
  m_contsurf_count = 0;

  //plastic
  malloc_t (pl_strain, double, m_elem_count * m_gp_count);
  malloc_t (sigma_y, double, m_elem_count * m_gp_count);  
  
  /// IMPLICIT THINGS
  //malloc_t(m_Kmat, Matrix, m_elem_count );  Written for asinglepointer
  m_Kmat = new Matrix*[m_elem_count];  // or use malloc_t macro if it's defined
  for (int e=0;e<m_elem_count;e++)
    m_Kmat[e] = new Matrix(m_nodxelem* m_dim,m_nodxelem* m_dim);

  
}

void Domain_d::AssignMaterial (Material_ *material_h) {
//    cudaMalloc((void**)&materials, 1 * sizeof(Material_ )); //
    printf("Assigning Material\n");
    malloc_t(materials, Material_,1);
    memcpy_t(materials, material_h, 1 * sizeof(Material_));	
}


  ///// ALREADY ALLOCATED
  //~ void Domain_d::setNode(const int &i, const double &_x, const double &_y, const double &_z){
    //~ if (i<m_node_count){
    //~ x[i*m_dim  ]=_x;
    //~ x[i*m_dim+1]=_y;
    //~ x[i*m_dim+2]=_z;
    //~ //return 1;
    //~ }
  //~ else{cout << "Node allocation error, node pos larger than node count."<<endl;}
        //~ //return 0;
  //~ }

};
