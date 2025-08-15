#include "ReMesher.h"
#include "Domain_d.h"

//#define DEBUG_MODE 0

#include "defs.h"

#include <iostream>
#include <limits>
#include <algorithm>


using namespace std;


namespace MetFEM{
  ReMesher::ReMesher(Domain_d *d)
  //:dim(d->m_dim)
  {

    m_dom = d;
    m_dom->m_dim = 3;

#ifdef REMESH_OMEGA_H 
    Omega_h::Library lib;
    Omega_h::Mesh mesh(&lib);


  #ifdef CUDA_BUILD
    // Allocate GPU memory
    double* d_node_data;
    int* d_connectivity_data;
    cudaMalloc((void**)&d_node_data, m_dom->m_node_count * 3 * sizeof(double));
    cudaMalloc((void**)&d_connectivity_data, num_elements * 4 * sizeof(int));

    // Copy from host to GPU
    //std::cout << "Old Mesh verts : "<<m_dom->m_node_count<<std::endl;
    cudaMemcpy(d_node_data, x, num_nodes * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_connectivity_data, h_connectivity_data, num_elements * 4 * sizeof(int), cudaMemcpyHostToDevice);

    // Create mesh using GPU data
    create_mesh(mesh, d_node_data, num_nodes, d_connectivity_data, num_elements);

    
    //vtk::Reader vtk_reader("piece_0.vtu");
    //classify_vertices(&mesh);
    //vtk_reader.read(mesh);

    // Free GPU memory
    cudaFree(d_node_data);
    cudaFree(d_connectivity_data);
  #else
    #ifdef DEBUG_MODE
    std::cout << "Creating mesh"<<std::endl;
    #endif
    // CPU case
    create_mesh(mesh, m_dom->x, m_dom->m_node_count, (int *)m_dom->m_elnod, m_dom->m_elem_count);
  #endif

      m_old_mesh = mesh; 
      m_mesh = mesh;


#endif //OMEGA_H

}

 
#include <array>
#include <iostream>

// Function to compute barycentric coordinates for a 3D tetrahedron
std::array<double, 4> barycentric_coordinates(const std::array<double, 3>& p,
                                              const std::array<double, 3>& p0,
                                              const std::array<double, 3>& p1,
                                              const std::array<double, 3>& p2,
                                              const std::array<double, 3>& p3) {
    // Compute volume of the tetrahedron
    std::array<double, 3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<double, 3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    std::array<double, 3> v2 = {p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};
    
    double detT = v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
                - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
                + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);

 if (std::abs(detT) < 1e-12) {  // Tolerance to avoid floating-point issues
        std::cerr << "Degenerate tetrahedron encountered!" << std::endl;
        return {-1.0, -1.0, -1.0, -1.0}; // Invalid coordinates
    }

    // Compute the barycentric coordinates
    std::array<double, 3> vp = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};

    double d0 = vp[0] * (v1[1] * v2[2] - v1[2] * v2[1])
              - vp[1] * (v1[0] * v2[2] - v1[2] * v2[0])
              + vp[2] * (v1[0] * v2[1] - v1[1] * v2[0]);

    double d1 = v0[0] * (vp[1] * v2[2] - vp[2] * v2[1])
              - v0[1] * (vp[0] * v2[2] - vp[2] * v2[0])
              + v0[2] * (vp[0] * v2[1] - vp[1] * v2[0]);

    double d2 = v0[0] * (v1[1] * vp[2] - v1[2] * vp[1])
              - v0[1] * (v1[0] * vp[2] - v1[2] * vp[0])
              + v0[2] * (v1[0] * vp[1] - v1[1] * vp[0]);

    double d3 = v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
              - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
              + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);

    double w0 = d0 / detT;
    double w1 = d1 / detT;
    double w2 = d2 / detT;
    double w3 = 1.0 - (w0 + w1 + w2);  // Enforce sum-to-one constraint

    return {w0, w1, w2, w3};
}

std::array<double, 4> stable_barycentric(const std::array<double, 3>& p,
                                        const std::array<double, 3>& p0,
                                        const std::array<double, 3>& p1,
                                        const std::array<double, 3>& p2,
                                        const std::array<double, 3>& p3) {
    const double EPS = 1.0e-9;
    const double NEG_EPS = -1.0e-7; // Slightly more tolerant for negatives
    
    // Compute vectors
    auto vec = [](const auto& a, const auto& b) {
        return std::array<double, 3>{b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    };
    
    std::array<double, 3> v0 = vec(p0, p1);
    std::array<double, 3> v1 = vec(p0, p2);
    std::array<double, 3> v2 = vec(p0, p3);
    std::array<double, 3> vp = vec(p0, p);

    // Compute determinant using more stable approach
    double det = v0[0]*(v1[1]*v2[2] - v1[2]*v2[1]) 
               - v0[1]*(v1[0]*v2[2] - v1[2]*v2[0]) 
               + v0[2]*(v1[0]*v2[1] - v1[1]*v2[0]);

    if (std::abs(det) < EPS) {
        return {NEG_EPS, NEG_EPS, NEG_EPS, NEG_EPS};
    }

    // Compute weights using the same basis for consistency
    double w0 = (vp[0]*(v1[1]*v2[2] - v1[2]*v2[1]) 
              - vp[1]*(v1[0]*v2[2] - v1[2]*v2[0]) 
              + vp[2]*(v1[0]*v2[1] - v1[1]*v2[0])) / det;

    double w1 = (v0[0]*(vp[1]*v2[2] - vp[2]*v2[1]) 
              - v0[1]*(vp[0]*v2[2] - vp[2]*v2[0]) 
              + v0[2]*(vp[0]*v2[1] - vp[1]*v2[0])) / det;

    double w2 = (v0[0]*(v1[1]*vp[2] - v1[2]*vp[1]) 
              - v0[1]*(v1[0]*vp[2] - v1[2]*vp[0]) 
              + v0[2]*(v1[0]*vp[1] - v1[1]*vp[0])) / det;

    double w3 = 1.0 - (w0 + w1 + w2);

    // Apply smoothing for near-zero weights
    if (std::abs(w0) < EPS) w0 = 0;
    if (std::abs(w1) < EPS) w1 = 0;
    if (std::abs(w2) < EPS) w2 = 0;
    if (std::abs(w3) < EPS) w3 = 0;

    return {w0, w1, w2, w3};
}

// Function to interpolate scalar values at the nodes of a tetrahedron
double interpolate_scalar(const std::array<double, 3>& p,
                          const std::array<double, 3>& p0, const std::array<double, 3>& p1, 
                          const std::array<double, 3>& p2, const std::array<double, 3>& p3,
                          double scalar0, double scalar1, double scalar2, double scalar3) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2, p3);
    return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2 + lambdas[3] * scalar3;
}

// Function to interpolate vector values at the nodes of a tetrahedron
std::array<double, 3> interpolate_vector(const std::array<double, 3>& p,
                                         const std::array<double, 3>& p0, const std::array<double, 3>& p1,
                                         const std::array<double, 3>& p2, const std::array<double, 3>& p3,
                                         std::array<double, 3> v0, std::array<double, 3> v1,
                                         std::array<double, 3> v2, std::array<double, 3> v3) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2, p3);
    return {
        lambdas[0] * v0[0] + lambdas[1] * v1[0] + lambdas[2] * v2[0] + lambdas[3] * v3[0],
        lambdas[0] * v0[1] + lambdas[1] * v1[1] + lambdas[2] * v2[1] + lambdas[3] * v3[1],
        lambdas[0] * v0[2] + lambdas[1] * v1[2] + lambdas[2] * v2[2] + lambdas[3] * v3[2]
    };
}

void setVec (double *vec, int size, int dim, double val = 0.0){

for (int v=0;v<size;v++)
    vec[v] = val;
  
}

void ReMesher::MapNodal(double *vfield, double *o_field){
  #ifndef REMESH_OMEGA_H
  MapNodalVectorRaw<3>        (vfield, o_field);
  #else              
  MapNodalVector<3> (m_mesh, vfield, o_field); 
  #endif
  
}
void ReMesher::MapNodalScalar(double *vfield, double *o_field){
  #ifndef REMESH_OMEGA_H
  MapNodalVectorRaw<1>        (vfield, o_field);
  #else              

  #endif
  
}

void ReMesher::MapElem(double *vfield, double *o_field, int field_dim){
  #ifndef REMESH_OMEGA_H ////MMG
  MapElemVectorRaw        (vfield, o_field,field_dim);
  #else              
  MapElemVector<3> (m_mesh, vfield, o_field,field_dim); 
 #endif
  
}

double tet_volume(double v0[3], double v1[3], double v2[3],  double v3[3]) {
    double a[3], b[3], c[3];
    for (int i = 0; i < 3; ++i) {
        a[i] = v1[i] - v0[i];
        b[i] = v2[i] - v0[i];
        c[i] = v3[i] - v0[i];
    }
    double volume = (a[0]*(b[1]*c[2] - b[2]*c[1])
                   - a[1]*(b[0]*c[2] - b[2]*c[0])
                   + a[2]*(b[0]*c[1] - b[1]*c[0])) / 6.0;
    return std::abs(volume);
}

//~ void LaplacianFilter(double* field, int* elnod, int elem_count, int nodxelem) {
  //~ std::vector<double> sum_field(m_node_count, 0.0);
  //~ std::vector<int> count(m_node_count, 0);

  //~ for (int e = 0; e < elem_count; e++) {
    //~ for (int a = 0; a < nodxelem; a++) {
      //~ int i = elnod[e * nodxelem + a];
      //~ sum_field[i] += field[i];
      //~ count[i]++;
    //~ }
  //~ }

  //~ for (int i = 0; i < m_node_count; i++) {
    //~ if (count[i] > 0)
      //~ field[i] = 0.5 * field[i] + 0.5 * (sum_field[i] / count[i]); // blend original + average
  //~ }
//~ }

void ReMesher::WriteDomain(){
  
  bool m_map_momentum = true;
  
  cout << "WRITING DOMAIN "<<m_node_count<<" NODES "<<m_elem_count<<"ELEMS"<<endl;  
  #ifdef REMESH_OMEGA_H
    m_node_count = m_mesh.nverts();
    m_elem_count = m_mesh.nelems();
  #endif

  //memcpy_t(m_->m_elnod, elnod_h, sizeof(int) * dom->m_elem_count * m_dom->m_nodxelem); 
  double *ufield  = new double [3*m_node_count];   
  //std::vector<double> ufield(3*m_node_count, 0.0);  
  double *vfield   = new double [3*m_node_count]; 
  double *afield   = new double [3*m_node_count]; 
  double *pafield   = new double [3*m_node_count];  //prev_a
  
  double *cforce   = new double [3*m_node_count];   
  
  double *esfield  = new double [m_elem_count]; 
  double *pfield   = new double [m_elem_count]; 
  double *sigfield = new double [6*m_elem_count]; 
  double *syfield  = new double [m_elem_count]; 
  double *psfield  = new double [m_elem_count];  
  double *str_rate = new double [6*m_elem_count]; 
  double *rot_rate = new double [6*m_elem_count];  
  double *tau      = new double [6*m_elem_count];   
  double *rho   = new double [m_elem_count];   
  double *vol_0 = new double [m_elem_count];    
  double *idetF   = new double [m_dom->m_elem_count];  //Inverse Deformation gradient
  
  double *Tfield  = new double [m_dom->m_node_count];
  
  double rho_0 = m_dom->rho_0[0];
  
  
  for (int i=0;i<3*m_node_count;i++){
    afield[i]=0.0;
    pafield[i]=0.0;
    vfield[i]=0.0;
    cforce[i] = 0.0;
  }
//  cout << "MAPPING"<<endl;
  ///// BEFORE REDIMENSION!
  //MapNodalVector<3>(m_mesh, ufield,  m_dom->u);
  cout <<"MAP NODAL"<<endl;
  //MapNodal(ufield,  m_dom->u);
  MapNodal(ufield,   m_dom->u); //new , old
  /////// IF MAP VEL DIRECTLY
  
  //if (m_dom->m_remesh_map_vel)
  ///// ATTENTION, MAP VELOCITY IS MANDATORY; IS BETTER TO MAP IT
  ///// THE PROBLEM IS THAT IS BENEFIT A LITTLE DAMPING
  if (!m_map_momentum)
    MapNodal(vfield,   m_dom->v); //DOES NOT CONS MOMENTUM
  
  if (m_dom->m_remesh_map_acc){
    MapNodal(afield,   m_dom->a);
    MapNodal(pafield,  m_dom->prev_a);
  }
  
  cout <<"DONE"<<endl;
  double *volumes=new double[m_elem_count];
  double vol = 0.0;
  for (int i=0;i<m_elem_count;i++){
      int n0 = m_elnod[4 * i];
      int n1 = m_elnod[4 * i + 1];
      int n2 = m_elnod[4 * i + 2];
      int n3 = m_elnod[4 * i + 3];
      double p0 []= {m_x[3*n0], m_x[3*n0+1], m_x[3*n0+2]};
      double p1 [] = {m_x[3*n1], m_x[3*n1+1], m_x[3*n1+2]};
      double p2 []= {m_x[3*n2], m_x[3*n2+1], m_x[3*n2+2]};
      double p3 [] = {m_x[3*n3], m_x[3*n3+1], m_x[3*n3+2]};
      volumes[i]=tet_volume(p0,p1,p2,p3);
      vol+=volumes[i];
  }
  cout << "New Volume: "<<vol<<endl;
  
  FindMapElemClosest();
  cout << "foundelem closrsrt"<<endl;
  MapElem(esfield,   m_dom->pl_strain);
  cout << "map pressure"<<endl;
  MapElem(pfield,    m_dom->p);

  // //cout << "mapping press"<<endl;
  // //HybridProjectionElemToElem<3>(m_mesh, pfield,    m_dom->p);

   MapElem(rho,    m_dom->rho);
   
  for (int e=0;e<m_dom->m_elem_count;e++)
    idetF[e] = m_dom->vol_0[e]/m_dom->vol[e];
  MapElem(vol_0,    idetF);
  
  
  // cout << "New mesh volume "<<total_volume<<endl;
  
   MapElem(sigfield,  m_dom->m_sigma         , 6);
   MapElem(str_rate,  m_dom->m_str_rate      , 6);
   MapElem(rot_rate,  m_dom->m_rot_rate      , 6);
   MapElem(tau,       m_dom->m_tau           , 6);
   
   MapElem(syfield,  m_dom->sigma_y  , 1);
   MapElem(psfield,  m_dom->pl_strain  , 1);

  cout << "Map done"<<endl;  
  // /////////////////////// BOUNDARY CONDITIONS
  int bccount[3];
  int    *bcx_nod,*bcy_nod,*bcz_nod;
  double *bcx_val,*bcy_val,*bcz_val;
  for (int d=0;d<3;d++) bccount[d]=m_dom->bc_count[d];
    
  bcx_nod =new int[bccount[0]];bcx_val =new double[bccount[0]];
  bcy_nod =new int[bccount[1]];bcy_val =new double[bccount[0]]; 
  bcz_nod =new int[bccount[2]];bcz_val =new double[bccount[0]];
  
  
    
  for (int b=0;b<bccount[0];b++){bcx_nod[b]=m_dom->bcx_nod[b];bcx_val[b]=m_dom->bcx_val[b];}
  for (int b=0;b<bccount[0];b++){bcy_nod[b]=m_dom->bcy_nod[b];bcy_val[b]=m_dom->bcy_val[b];}
  for (int b=0;b<bccount[0];b++){bcz_nod[b]=m_dom->bcz_nod[b];bcz_val[b]=m_dom->bcz_val[b];}
  
  /// BEFORE DELETE, SAVE MASS FOR CONSERVATION
  if (m_map_momentum){
  double old_p_field[m_dom->m_dim*m_dom->m_node_count];
  for (int i=0;i<m_dom->m_node_count;i++)
    for (int d=0;d<m_dom->m_dim; d++)
      old_p_field[m_dom->m_dim*i+d] = m_dom->v[m_dom->m_dim*i+d]*m_dom->m_mdiag[i];
  
  ////IF MAP MOMENTUM
  MapNodal(vfield,   old_p_field); 
  }
  
  ////BEFORE REWRITE
  //// WRITE
  m_dom->Free();
  
   #ifdef REMESH_OMEGA_H
    m_dom->m_node_count = m_mesh.nverts();
    m_dom->m_elem_count = m_mesh.nelems();
  #else
    
    m_dom->m_node_count = m_node_count;
    m_dom->m_elem_count = m_elem_count;    
  #endif
  
  cout << "creating domain"<<endl;
  m_dom->SetDimension(m_dom->m_node_count,m_dom->m_elem_count);	 //AFTER CREATING DOMAIN

  malloc_t(m_dom->m_elnod, unsigned int,m_dom->m_elem_count * m_dom->m_nodxelem);

  
  cout << "SETTING DENSITY TO "<<rho_0<<endl;
  m_dom->setDensity(rho_0);


  // // ATENTION:
  // //Volumes could also be calculated with Jacobian
  for (int e=0;e<m_dom->m_elem_count;e++){
    m_dom->vol_0[e] = vol_0[e]*volumes[e];
    m_dom->vol  [e] = volumes[e];
  }
  
  ///// TO MOMENTUM CONSERVATION
    

  // //cout << "COPYING "<<m_dom->m_elem_count * m_dom->m_nodxelem<< " element nodes "<<endl;
  memcpy_t(m_dom->u,        ufield, sizeof(double) * m_dom->m_node_count * 3);  
  memcpy_t(m_dom->v,        vfield, sizeof(double) * m_dom->m_node_count * 3);
  memcpy_t(m_dom->m_vprev,  vfield, sizeof(double) * m_dom->m_node_count * 3);  
  
  memcpy_t(m_dom->a,        afield, sizeof(double) * m_dom->m_node_count * 3);   
  memcpy_t(m_dom->prev_a,  pafield, sizeof(double) * m_dom->m_node_count * 3);   
  memcpy_t(m_dom->contforce,cforce, sizeof(double) * m_dom->m_node_count * 3);

  memcpy_t(m_dom->pl_strain,      esfield,  sizeof(double) * m_dom->m_elem_count ); 
  memcpy_t(m_dom->pl_strain_prev, esfield,  sizeof(double) * m_dom->m_elem_count ); 
  
  memcpy_t(m_dom->m_sigma  ,      sigfield,   sizeof(double) * m_dom->m_elem_count *6); 
  memcpy_t(m_dom->m_sigma_prev,   sigfield,   sizeof(double) * m_dom->m_elem_count *6); //TO BLEND
  memcpy_t(m_dom->m_str_rate,     str_rate,   sizeof(double) * m_dom->m_elem_count *6); 
  memcpy_t(m_dom->m_str_rate_prev,str_rate,   sizeof(double) * m_dom->m_elem_count *6); 
  
  memcpy_t(m_dom->m_rot_rate,     rot_rate,   sizeof(double) * m_dom->m_elem_count *6); 
  memcpy_t(m_dom->m_tau,          tau,        sizeof(double) * m_dom->m_elem_count *6); 
  memcpy_t(m_dom->m_tau_prev,     tau,   sizeof(double) * m_dom->m_elem_count *6); //TO BLEND
    
  memcpy_t(m_dom->sigma_y,   syfield,  sizeof(double) * m_dom->m_elem_count ); 
  memcpy_t(m_dom->pl_strain, psfield,  sizeof(double) * m_dom->m_elem_count ); 

  memcpy_t(m_dom->p,          pfield,  sizeof(double) * m_dom->m_elem_count ); 
    
  memcpy_t(m_dom->rho,          rho,  sizeof(double) * m_dom->m_elem_count ); 
  
  if (m_dom->m_thermal){
    memcpy_t(m_dom->T,        Tfield, sizeof(double) * m_dom->m_node_count);     
    
  }

   
  m_dom->AssignMatAddress();
  const Material_ *matt  = &m_dom->materials[0];
  cout << "G "<<matt->Elastic().G()<<endl;
  
  ////// AFTER ALL COPIES //////
  cout << "copying"<<endl;
  memcpy_t(m_dom->x,      m_x, 3*sizeof(double) * m_dom->m_node_count);       
  memcpy_t(m_dom->m_elnod,  m_elnod, 4*sizeof(int) * m_dom->m_elem_count);  
    cout << "Setting conn"<<endl;
  m_dom->setNodElem(m_elnod);     

  //Now calculate new mass, after mapping volume
  cout << "Calculating masses"<<endl;
  m_dom->CalcNodalVol();
  m_dom->CalcNodalMassFromVol();
  //if (m_dom->m_remesh_map_vel)
 cout << "recovering velocities"<<endl;
  if (m_map_momentum){
    for (int i=0;i<m_node_count;i++){ //Or already dom_d->m_node_count since domain changed
      for (int d=0;d<m_dom->m_dim;d++)
        vfield[m_dom->m_dim*i+d] /=m_dom->m_mdiag[i];
    }
    memcpy_t(m_dom->m_vprev,  vfield, sizeof(double) * m_dom->m_node_count * 3); 
    memcpy_t(m_dom->v,        vfield, sizeof(double) * m_dom->m_node_count * 3); 
  }
  
  //// RECALCULATED FROM MOMENTUM
  ///////if (m_dom->m_remesh_map_vel)
  /////memcpy_t(m_dom->v,        vfield, sizeof(double) * m_dom->m_node_count * 3); 
          
    
  
  ReMapBCs(bcx_nod,bcx_val,m_dom->bcx_nod, m_dom->bcx_val, bccount[0]);
  ReMapBCs(bcy_nod,bcy_val,m_dom->bcy_nod, m_dom->bcy_val, bccount[1]);
  ReMapBCs(bcz_nod,bcz_val,m_dom->bcz_nod, m_dom->bcz_val, bccount[2]);

  
  cout << "deleting "<<endl;
  delete [] ufield, vfield, afield, pafield;
  delete [] esfield,pfield,sigfield, syfield, psfield;
  delete [] str_rate,rot_rate, tau,rho,vol_0,idetF;
  delete [] bcx_nod,bcy_nod,bcz_nod,bcx_val,bcy_val,bcz_val;
  delete [] Tfield;
  cout << "MESH CHANGED"<<endl;

  //AFTER MAP
  //THIS CRASHES
  //free_t(m_closest_elem);
  

}


/////WITH THE RAW ELEM AND CONNECT
//args: NEW (m_node_count), OLD(m_dom->m_elem_count)
template <int dim>
void ReMesher::MapNodalVectorRaw(double *vfield, double *o_field) {
    // Loop over the target nodes in the new mesh

    cout << "MAP NODAL VECTOR RAW (MMMG)"<<endl;
    for (int vert=0;vert<m_node_count;vert++){
      
        for (int d=0;d<3;d++) vfield [dim*vert+d] = 0.0;

   //auto coords = mesh.coords();

    
    //auto f = OMEGA_H_LAMBDA(LO vert) {
      //auto x = get_vector<3>(coords, vert);
      double x[3];
      for (int d=0;d<3;d++) x[d]=m_x[3*vert+d];
      
      bool found_samenode = false;
      double tol = 1.0e-6;
      //SEARCH OVERALL NEW MESH NODES IF NOT A NEW NODE NEAR THE OLD NODE  
      for (int v = 0; v < m_dom->m_node_count; v++){
        //If new node dist <tol, map new node = old node
        std::array<double, 3> p0 = {m_dom->x[3 * v], m_dom->x[3 * v + 1], m_dom->x[3 * v + 2]};
        double distance = 0.0;
        for (int i = 0; i < 3; ++i) {
            distance += (x[i] - p0[i]) * (x[i] - p0[i]);
        }
        if (distance<tol){
          found_samenode = true;
          //cout << "FOUND " << vert << " SAME NODE "<<endl;
          if (vert == 0)
            std::cout << "FOUND new node "<<v << " For node "<<vert<<std::endl;
          
          for (int d=0;d<dim;d++) vfield[dim*vert+d] = o_field[dim*v+d];
        }                
      }//node
      
      if (!found_samenode){
    //for (int n = 0; n < mesh.nverts(); n++) {
        bool found = false;  // Flag to indicate whether the node is inside an element in the old mesh
        //std::cout << mesh.coords()[n]<<std::endl;
        
        // Get coordinates for the node in the new mesh
        std::array<double, 3> target_node = {x[0], x[1], x[2]}; // Now using 3D coordinates
        
        // Loop over the elements in the old mesh (using *elnod to access connectivity and *node for coordinates)
        for (int i = 0; i < m_dom->m_elem_count; i++) {
            // Connectivity for the tetrahedral element (assumed to have 4 nodes per element in the old mesh)
            int n0 = m_dom->m_elnod[4*i];   // Node 0 in the element
            int n1 = m_dom->m_elnod[4*i+1]; // Node 1 in the element
            int n2 = m_dom->m_elnod[4*i+2]; // Node 2 in the element
            int n3 = m_dom->m_elnod[4*i+3]; // Node 3 in the element

            std::array<double, 3> p0 = {m_dom->x[3*n0], m_dom->x[3*n0+1], m_dom->x[3*n0+2]};
            std::array<double, 3> p1 = {m_dom->x[3*n1], m_dom->x[3*n1+1], m_dom->x[3*n1+2]};
            std::array<double, 3> p2 = {m_dom->x[3*n2], m_dom->x[3*n2+1], m_dom->x[3*n2+2]};
            std::array<double, 3> p3 = {m_dom->x[3*n3], m_dom->x[3*n3+1], m_dom->x[3*n3+2]};

            std::array<double, 4> lambdas = stable_barycentric(target_node, p0, p1, p2, p3);

            if (lambdas[0] >= -1.0e-8 && lambdas[1] >= -1.0e-8 && lambdas[2] >= -1.0e-8 && lambdas[3] >= -1.0e-8) { 

                //double scalar[4];
                //for (int n=0;n<4;n++) scalar[n] = m_dom->pl_strain[i];

                //double interpolated_scalar = interpolate_scalar(target_node, p0, p1, p2, p3, scalar[0], scalar[1], scalar[2], scalar[3]);


                // Interpolate vector values for displacement (if needed)
                std::array<double, 3> disp[4];
                for (int n=0;n<4;n++)
                  for (int d=0;d<dim;d++)
                    disp[n][d] = o_field[dim*m_dom->m_elnod[4*i+n]+d];
                
                //cout << "Interp disp"<<endl;
                std::array<double, 3> interpolated_disp = interpolate_vector(target_node, p0, p1, p2, p3, disp[0], disp[1], disp[2],disp[3]);
                for (int d=0;d<dim;d++) vfield[dim*vert+d] = interpolated_disp[d];
                // Optionally, interpolate other scalar/vector fields for the new mesh node here
                if (vert == 0)  {
                std::cout << "FOUND ELEMENT "<<i << " For node "<<vert<<std::endl;
                //std::cout << "Node " << vert << " is inside element " << i << " of the old mesh." << std::endl;
                  //std::cout << "Interpolated scalar: " << interpolated_scalar << std::endl;
                  std::cout << "Interpolated displacement: (" << interpolated_disp[0] << ", " << interpolated_disp[1] << ", " << interpolated_disp[2] << ")\n";
                }
                found = true;
                break;  // Exit the element loop once the element is found
            }//lambdas
          }//elem
          if (!found) {
              std::cout << "Node " << vert << " is not inside any element of the old mesh." << std::endl;
          }
        }//found same node

      //n++;
  
    }

}//MAP


int ReMesher::find_closest_node(const double x[3]) {
    int closest = -1;
    double min_dist_sq = 1.0e20;

    for (int i = 0; i < m_dom->m_node_count; ++i) {
        double dx = x[0] - m_dom->x[3*i];
        double dy = x[1] - m_dom->x[3*i + 1];
        double dz = x[2] - m_dom->x[3*i + 2];
        double dist_sq = dx*dx + dy*dy + dz*dz;

        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            closest = i;
        }
    }

    return closest;
}

// void ReMesher::MapNodalVectorRaw(double *vfield, double *o_field) {
    // const double EPS = 1.0e-4;
    // const double SAME_NODE_TOL = 1.0e-6;

    // // double avg_element_size = ...; // estimate from mesh, e.g., average edge length or cube root of average element volume
    // // double tolerance = avg_element_size * 1e-3; // e.g. 0.001 times element size

    // int fallback_node = 16;
    // double xx[3] = {m_x[3*fallback_node], m_x[3*fallback_node+1], m_x[3*fallback_node+2]};
    // bool found = false;

    // for (int i = 0; i < m_dom->m_elem_count; ++i) {
        // int n0 = m_dom->m_elnod[4*i];
        // int n1 = m_dom->m_elnod[4*i+1];
        // int n2 = m_dom->m_elnod[4*i+2];
        // int n3 = m_dom->m_elnod[4*i+3];

        // std::array<double, 3> p0 = {m_dom->x[3*n0], m_dom->x[3*n0+1], m_dom->x[3*n0+2]};
        // std::array<double, 3> p1 = {m_dom->x[3*n1], m_dom->x[3*n1+1], m_dom->x[3*n1+2]};
        // std::array<double, 3> p2 = {m_dom->x[3*n2], m_dom->x[3*n2+1], m_dom->x[3*n2+2]};
        // std::array<double, 3> p3 = {m_dom->x[3*n3], m_dom->x[3*n3+1], m_dom->x[3*n3+2]};

        // auto lambdas = stable_barycentric({xx[0], xx[1], xx[2]}, p0, p1, p2, p3);
        // if (lambdas[0] >= -1e-5 && lambdas[1] >= -1e-5 &&
            // lambdas[2] >= -1e-5 && lambdas[3] >= -1e-5) {
            // std::cout << "✅ Node " << fallback_node << " lies inside element " << i << std::endl;
            // found = true;
            // break;
        // }
    // }
    // if (!found){
      // cout <<"Node 16 not found inside elements!!"<<endl;
      
    // }
    
    // // Pre-compute element centers for faster searching
    // std::vector<std::array< double, 3> > elem_centers(m_dom->m_elem_count);
    // for (int i = 0; i < m_dom->m_elem_count; i++) {
        // elem_centers[i] = {0,0,0};
        // for (int n = 0; n < 4; n++) {
            // int node_idx = m_dom->m_elnod[4*i + n];
            // elem_centers[i][0] += m_dom->x[3*node_idx];
            // elem_centers[i][1] += m_dom->x[3*node_idx + 1];
            // elem_centers[i][2] += m_dom->x[3*node_idx + 2];
        // }
        // for (int d = 0; d < 3; d++) elem_centers[i][d] /= 4.0;
    // }

    // for (int vert = 0; vert < m_dom->m_node_count; vert++) {
        // double x[3] = {m_x[3*vert], m_x[3*vert+1], m_x[3*vert+2]};
        
        // // 1. Check for exact node matches first
        // bool found_samenode = false;
        // for (int v = 0; v < m_dom->m_node_count; v++) {
            // double dist_sq = 0.0;
            // for (int d = 0; d < 3; d++) {
                // double diff = x[d] - m_dom->x[3*v+d];
                // dist_sq += diff*diff;
            // }
            // if (dist_sq < SAME_NODE_TOL*SAME_NODE_TOL) {
                // for (int d = 0; d < 3; d++) 
                    // vfield[3*vert+d] = o_field[3*v+d];
                // found_samenode = true;
                // break;
            // }
        // }
        // if (found_samenode) continue;

        // // 2. Find closest elements (for smoother transitions)
        // const int NUM_CANDIDATES = 20;
        // std::vector<std::pair<double, int>> candidate_elements;
        
        // for (int i = 0; i < m_dom->m_elem_count; i++) {
            // double dist_sq = 0.0;
            // for (int d = 0; d < 3; d++) {
                // double diff = x[d] - elem_centers[i][d];
                // dist_sq += diff*diff;
            // }
            // candidate_elements.emplace_back(dist_sq, i);
        // }
        
        // // Sort by distance and take top candidates
        // std::sort(candidate_elements.begin(), candidate_elements.end());
        // candidate_elements.resize(std::min(NUM_CANDIDATES, (int)candidate_elements.size()));

        // // 3. Try each candidate element for valid interpolation
        // bool found = false;
        // int found_el;
        // std::array<double, 3> interpolated_disp = {0,0,0};
        
        // for (const auto& [dist, i] : candidate_elements) {
            // int n0 = m_dom->m_elnod[4*i];
            // int n1 = m_dom->m_elnod[4*i+1];
            // int n2 = m_dom->m_elnod[4*i+2];
            // int n3 = m_dom->m_elnod[4*i+3];

            // std::array<double, 3> p0 = {m_dom->x[3*n0], m_dom->x[3*n0+1], m_dom->x[3*n0+2]};
            // std::array<double, 3> p1 = {m_dom->x[3*n1], m_dom->x[3*n1+1], m_dom->x[3*n1+2]};
            // std::array<double, 3> p2 = {m_dom->x[3*n2], m_dom->x[3*n2+1], m_dom->x[3*n2+2]};
            // std::array<double, 3> p3 = {m_dom->x[3*n3], m_dom->x[3*n3+1], m_dom->x[3*n3+2]};
        
            // auto lambdas = stable_barycentric({x[0],x[1],x[2]}, p0, p1, p2, p3);

        
            // if (lambdas[0] >= -EPS && lambdas[1] >= -EPS && 
                // lambdas[2] >= -EPS && lambdas[3] >= -EPS) {
                
                // std::array<double, 3> disp[4];
                // for (int n = 0; n < 4; n++) {
                    // int node_idx = m_dom->m_elnod[4*i+n];
                    // disp[n] = {o_field[3*node_idx], o_field[3*node_idx+1], o_field[3*node_idx+2]};
                // }
                
                // // Weighted average that preserves constant Z-values
                // interpolated_disp = {0,0,0};
                // for (int n = 0; n < 4; n++) {
                    // double w = lambdas[n];
                    // if (w < 0) w = 0; // Clip negative weights
                    // interpolated_disp[0] += w * disp[n][0];
                    // interpolated_disp[1] += w * disp[n][1];
                    // interpolated_disp[2] += w * disp[n][2];
                // }
                
                // found = true;
                // found_el = i;
                // break;
            // }
        // }

        // if (found) {
            // for (int d = 0; d < 3; d++) 
                // vfield[3*vert+d] = interpolated_disp[d];
            // //cout << "Node "<<vert <<"inside element "<<found_el<<endl;
        // } else {
            // // Fallback: Use nearest node value
            // int closest_node = find_closest_node(x);
            // for (int d = 0; d < 3; d++)
                // vfield[3*vert+d] = o_field[3*closest_node+d];
            // std::cout << "Warning: Using fallback for node " << vert << std::endl;
        // }
    // }
// }

// THIS NOT USES MESH AS OUTPUT; USES 
////DO THIS AFTER SERACH CLOSEST
void ReMesher::MapElemVectorRaw(double *vfield, double *o_field, int field_dim) {
  
  for (int elem=0;elem<m_elem_count;elem++){ ///LOOP THROUGH NEW MESH  CELLS
      if(m_closest_elem[elem]>-1)  
        for (int d=0;d<field_dim;d++) 
          vfield[elem*field_dim+d] = o_field[m_closest_elem[elem]*field_dim+d];  
      else 
        for (int d=0;d<field_dim;d++) vfield[elem*field_dim+d]=0.0;  
  }
}

/////////////////////////////
////// OLD FUNCTION /////////
//// FIND CLOSEST ////

//~ void ReMesher::FindMapElemClosest() {


    //~ for (int elem=0;elem<m_elem_count;elem++){ ///LOOP THROUGH NEW MESH  CELLS
    //~ std::array<double, 3> barycenter = {0.0, 0.0, 0.0};

    //~ std::array<double, 3> barycenter_old_clos = {0.0, 0.0, 0.0};
    //~ bool found = false;
        //~ // Calculate barycenter of the current new element
        //~ for (int en = 0; en < 4; en++) {
            //~ //auto v = elems2verts.ab2b[elem * 4 + en];
            //~ //auto x = get_vector<3>(coords, v);
            //~ int v = m_elnod[elem * 4 + en];

            //~ double x[3];
            //~ for (int d=0;d<3;d++)x[d]=m_x[3*v+d]; //X: NEW MESH NODE COORDS

            //~ barycenter[0] += x[0];
            //~ barycenter[1] += x[1];
            //~ barycenter[2] += x[2];
        //~ }
        //~ barycenter[0] /= 4.0;
        //~ barycenter[1] /= 4.0;
        //~ barycenter[2] /= 4.0;

        //~ // Search for the closest old element by distance
        //~ double min_distance = std::numeric_limits<double>::max();
        //~ int closest_elem = -1;

        //~ for (int i = 0; i < m_dom->m_elem_count; i++) {
            //~ int n0 = m_dom->m_elnod[4 * i];
            //~ int n1 = m_dom->m_elnod[4 * i + 1];
            //~ int n2 = m_dom->m_elnod[4 * i + 2];
            //~ int n3 = m_dom->m_elnod[4 * i + 3];

            //~ std::array<double, 3> p0 = {m_dom->x[3 * n0], m_dom->x[3 * n0 + 1], m_dom->x[3 * n0 + 2]};
            //~ std::array<double, 3> p1 = {m_dom->x[3 * n1], m_dom->x[3 * n1 + 1], m_dom->x[3 * n1 + 2]};
            //~ std::array<double, 3> p2 = {m_dom->x[3 * n2], m_dom->x[3 * n2 + 1], m_dom->x[3 * n2 + 2]};
            //~ std::array<double, 3> p3 = {m_dom->x[3 * n3], m_dom->x[3 * n3 + 1], m_dom->x[3 * n3 + 2]};

            //~ // Calculate the barycenter of the old element
            //~ std::array<double, 3> old_barycenter = {
                //~ (p0[0] + p1[0] + p2[0] + p3[0]) / 4.0,
                //~ (p0[1] + p1[1] + p2[1] + p3[1]) / 4.0,
                //~ (p0[2] + p1[2] + p2[2] + p3[2]) / 4.0
            //~ };

            //~ double distance = 
            //~ //std::sqrt(
                //~ std::pow(barycenter[0] - old_barycenter[0], 2) +
                //~ std::pow(barycenter[1] - old_barycenter[1], 2) +
                //~ std::pow(barycenter[2] - old_barycenter[2], 2)
            //~ ;
            //~ //);

            //~ if (distance < min_distance) {
                //~ min_distance = distance;
                //~ m_closest_elem[elem] = i;
                //~ found = true;
                //~ barycenter_old_clos = old_barycenter;
            //~ }
        //~ }//elem

        //~ if (found) {
            //~ //std::cout << "Mapped element " << elem << " to old element " << closest_elem << std::endl;
        //~ } else {
            //~ std::cout << "ERROR: No matching element found for element " << elem << std::endl;
        //~ }

    //~ }//elem
//~ }

/////////////////////NEW 

void ReMesher::FindMapElemClosest() {
    for (int elem = 0; elem < m_elem_count; elem++) {
        std::array<double, 3> barycenter = {0.0, 0.0, 0.0};
        
        //////Calculate barycenter of current new element
        for (int en = 0; en < 4; en++) {
            int v = m_elnod[elem * 4 + en];
            barycenter[0] += m_x[3 * v];
            barycenter[1] += m_x[3 * v + 1];
            barycenter[2] += m_x[3 * v + 2];
        }
        barycenter[0] /= 4.0;
        barycenter[1] /= 4.0;
        barycenter[2] /= 4.0;

        bool found_containing = false;
        int closest_elem = -1;
        double min_distance = std::numeric_limits<double>::max();
        std::array<double, 3> closest_barycenter;

        ////////First pass: Try to find an element that contains the barycenter
        for (int i = 0; i < m_dom->m_elem_count; i++) {
            int n0 = m_dom->m_elnod[4 * i];
            int n1 = m_dom->m_elnod[4 * i + 1];
            int n2 = m_dom->m_elnod[4 * i + 2];
            int n3 = m_dom->m_elnod[4 * i + 3];

            std::array<double, 3> p0 = {m_dom->x[3 * n0], m_dom->x[3 * n0 + 1], m_dom->x[3 * n0 + 2]};
            std::array<double, 3> p1 = {m_dom->x[3 * n1], m_dom->x[3 * n1 + 1], m_dom->x[3 * n1 + 2]};
            std::array<double, 3> p2 = {m_dom->x[3 * n2], m_dom->x[3 * n2 + 1], m_dom->x[3 * n2 + 2]};
            std::array<double, 3> p3 = {m_dom->x[3 * n3], m_dom->x[3 * n3 + 1], m_dom->x[3 * n3 + 2]};

            ////////Reuse your existing barycentric coordinate function
            std::array<double, 4> lambdas = stable_barycentric(barycenter, p0, p1, p2, p3);

            if (lambdas[0] >= -1.0e-8 && lambdas[1] >= -1.0e-8 && 
                lambdas[2] >= -1.0e-8 && lambdas[3] >= -1.0e-8) {
                m_closest_elem[elem] = i;
                found_containing = true;
                
                ////////Diagnostic output
                if (elem < 5) { ////////Print first few elements for verification
                    std::cout << "Element " << elem << " barycenter INSIDE old element " << i 
                              << " with barycentric coords: (" 
                              << lambdas[0] << ", " << lambdas[1] << ", " 
                              << lambdas[2] << ", " << lambdas[3] << ")\n";
                }
                break;
            }

            /////////////////If not inside, compute distance for fallback
            std::array<double, 3> old_barycenter = {
                (p0[0] + p1[0] + p2[0] + p3[0]) / 4.0,
                (p0[1] + p1[1] + p2[1] + p3[1]) / 4.0,
                (p0[2] + p1[2] + p2[2] + p3[2]) / 4.0
            };

            double distance = std::pow(barycenter[0] - old_barycenter[0], 2) +
                             std::pow(barycenter[1] - old_barycenter[1], 2) +
                             std::pow(barycenter[2] - old_barycenter[2], 2);

            if (distance < min_distance) {
                min_distance = distance;
                closest_elem = i;
                closest_barycenter = old_barycenter;
            }
        }

        ///////////////Fallback to closest element if no containing element found
        if (!found_containing) {
            m_closest_elem[elem] = closest_elem;
            
            ////////////////Diagnostic output
            if (elem < 5) {
                std::cout << "Element " << elem << " barycenter NOT INSIDE any element. "
                          << "Using closest element " << closest_elem 
                          << " with distance " << std::sqrt(min_distance) << "\n";
                std::cout << "New barycenter: (" << barycenter[0] << ", " 
                          << barycenter[1] << ", " << barycenter[2] << ")\n";
                std::cout << "Old barycenter: (" << closest_barycenter[0] << ", " 
                          << closest_barycenter[1] << ", " << closest_barycenter[2] << ")\n";
            }
        }
    }
}


// 1. Using the keep_classification Option
// The simplest way to preserve certain entities is during adaptation setup:

// cpp
// Copy
// Omega_h::AdaptOpts opts(&mesh);
// opts.keep_classification = true; // Preserves classification tags
// opts.verbosity = Omega_h::EXTRA_STATS;
// Omega_h::adapt(&mesh, opts);
// 2. Explicit Mapping with TransferPair
// For precise tracking of unchanged entities:

// cpp
// Copy
// // Before adaptation
// Omega_h::LOs old_verts = mesh.ask_verts_of_dim(mesh.dim());
// Omega_h::LOs old_elems = mesh.ask_elements_of_dim(mesh.dim());

// // Perform adaptation
// Omega_h::adapt(&mesh, opts);

// // Get mapping after adaptation
// auto verts_transfer = mesh.get_transfer_map(Omega_h::VERT);
// auto elems_transfer = mesh.get_transfer_map(Omega_h::REGION);

// // Identify unchanged vertices
// Omega_h::Write<Omega_h::LO> unchanged_verts(old_verts.size(), 0);
// Omega_h::parallel_for(old_verts.size(), OMEGA_H_LAMBDA(Omega_h::LO v) {
    // if (verts_transfer[v] != Omega_h::LO(-1)) {
        // unchanged_verts[v] = 1; // Mark as unchanged
    // }
// });

void ReMesher::ReMapBCs(int  *old_bc_nod,
                      double *old_bc_val,

                    int  *new_bc_nod,
                    double *new_bc_val,
                    int bc_count) {
    
  // //cout << "MAPPING BCs"<<endl;

  // //cout << "coords"<<endl;
  // auto new_coords = m_mesh.coords();
  // //cout << "OLDCOORDS"<<endl;
  // auto old_coords = m_old_mesh.coords();
  // //int dim = m_old_mesh.dim();
  
  // int new_count = 0;
  // for (std::size_t i = 0; i < bc_count; ++i) {
    // //cout << "vert "<<i<<endl;
    // I64 old_id = old_bc_nod[i];
    // Real val = old_bc_val[i];

    // Vector<3> p_old;
    // for (int d = 0; d < 3; ++d) {
      // p_old[d] = old_coords[old_id * 3 + d];
    // }

    // Real min_dist2 = std::numeric_limits<Real>::max();
    // I64 closest_id = -1;

    // for (I64 new_id = 0; new_id < m_mesh.nverts(); ++new_id) {
      // Vector<3> p_new;
      // for (int d = 0; d < 3; ++d) {
        // p_new[d] = new_coords[new_id * 3 + d];
      // }

      // Real dist2 = norm_squared(p_new - p_old);
      // if (dist2 < min_dist2) {
        // min_dist2 = dist2;
        // closest_id = new_id;
      // }
    // }

    // if (closest_id >= 0) {
      // if (closest_id != i)
        // //cout << "Different Node id found ,old "<<i<<", New "<< closest_id <<endl;
      // new_bc_nod[i]=closest_id;
      // new_bc_val[i]=val;
      // //cout << "val "<<val<<endl;
    // } else {
      // std::cerr << "Warning: Could not find nearest node for BC node " << old_id << std::endl;
    // }
  // }//bc_count
}




void ReMesher::ReMapBCsByFace(int* old_bc_nod,
                        double* old_bc_val,
                        int* new_bc_nod,
                        double* new_bc_val,
                        int bc_count) {
  
  // auto new_coords = m_mesh.coords();
  // auto old_coords = m_old_mesh.coords();
  // auto elems2verts = m_mesh.ask_down(3, 0);
  // int nelems = m_mesh.nelems();

  // // Step 1: Build surface triangle list (triangles with only one adjacent tet)
  // std::map<std::array<int, 3>, int> face_count;
  // for (int e = 0; e < nelems; ++e) {
    // int n0 = elems2verts.ab[4 * e + i]
    // int n1 = elems2verts.ab[4 * e + 1];
    // int n2 = elems2verts.ab[4 * e + 2];
    // int n3 = elems2verts.ab[4 * e + 3];
    
    // int faces[4][3] = {
      // {n0, n1, n2}, {n0, n1, n3},
      // {n0, n2, n3}, {n1, n2, n3}
    // };
    
    // for (int f = 0; f < 4; ++f) {
      // std::array<int, 3> tri = {faces[f][0], faces[f][1], faces[f][2]};
      // std::sort(tri.begin(), tri.end());
      // face_count[tri]++;
    // }
  // }

  // std::vector<std::array<int, 3>> boundary_faces;
  // for (auto& [face, count] : face_count) {
    // if (count == 1) boundary_faces.push_back(face);
  // }

  // // Step 2: Map each BC point onto boundary faces
  // for (int i = 0; i < bc_count; ++i) {
    // int old_id = old_bc_nod[i];
    // Real val = old_bc_val[i];

    // Vector<3> p_old;
    // for (int d = 0; d < 3; ++d) {
      // p_old[d] = old_coords[old_id * 3 + d];
    // }

    // Real min_dist2 = std::numeric_limits<Real>::max();
    // std::array<int, 3> best_face = {-1, -1, -1};
    // Vector<3> best_bary = {0, 0, 0};

    // for (const auto& face : boundary_faces) {
      // Vector<3> a, b, c;
      // for (int d = 0; d < 3; ++d) {
        // a[d] = new_coords[face[0] * 3 + d];
        // b[d] = new_coords[face[1] * 3 + d];
        // c[d] = new_coords[face[2] * 3 + d];
      // }

      // // Compute normal and barycentric coordinates
      // Vector<3> v0 = b - a;
      // Vector<3> v1 = c - a;
      // Vector<3> v2 = p_old - a;

      // Real d00 = dot(v0, v0);
      // Real d01 = dot(v0, v1);
      // Real d11 = dot(v1, v1);
      // Real d20 = dot(v2, v0);
      // Real d21 = dot(v2, v1);
      // Real denom = d00 * d11 - d01 * d01;

      // if (std::abs(denom) < 1e-12) continue;

      // Real v = (d11 * d20 - d01 * d21) / denom;
      // Real w = (d00 * d21 - d01 * d20) / denom;
      // Real u = 1.0 - v - w;

      // if (u >= -1e-4 && v >= -1e-4 && w >= -1e-4) {
        // Vector<3> proj = u * a + v * b + w * c;
        // Real dist2 = norm_squared(p_old - proj);
        // if (dist2 < min_dist2) {
          // min_dist2 = dist2;
          // best_face = face;
          // best_bary = {u, v, w};
        // }
      // }
    // }

    // if (best_face[0] != -1) {
      // // Here: assign value to node with max barycentric weight
      // int max_idx = 0;
      // if (best_bary[1] > best_bary[max_idx]) max_idx = 1;
      // if (best_bary[2] > best_bary[max_idx]) max_idx = 2;

      // int best_node = best_face[max_idx];
      // new_bc_nod[i] = best_node;
      // new_bc_val[i] = val;
    // } else {
      // std::cerr << "Warning: No triangle match found for BC node " << old_id << std::endl;
    // }
  // }
  
  
}
  
};
