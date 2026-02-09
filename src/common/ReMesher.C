/*************************************************************************/
/*  ReMesher.C                                                   */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                */
/*  ('https://www.opensourcemech.com',)                                    */
/*                                                                       */
/*  Copyright (c) 2023-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. */ 
/*  See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/



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

    m_params.nodal_tol    = 1.0e-3; //relative to elem length
    m_params.lambda_tol   = 1.0e-4; //relative to elem length

}

 
#include <array>
#include <iostream>

double tri_area(const double p0[2], const double p1[2], const double p2[2]) {
    return 0.5 * std::abs(
        (p1[0] - p0[0]) * (p2[1] - p0[1]) -
        (p1[1] - p0[1]) * (p2[0] - p0[0])
    );
}

double quad_area(const double p0[2],
                 const double p1[2],
                 const double p2[2],
                 const double p3[2]) {
    // Split: (0,1,2) + (0,2,3)
    return tri_area(p0, p1, p2) + tri_area(p0, p2, p3);
}


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

//~ template <int dim>
//~ std::array<double, dim>
//~ interpolate_vector(
    //~ const std::array<double, 3>& p,
    //~ const std::array<double, 3>& p0,
    //~ const std::array<double, 3>& p1,
    //~ const std::array<double, 3>& p2,
    //~ const std::array<double, 3>& p3,
    //~ const std::array<double, dim>& v0,
    //~ const std::array<double, dim>& v1,
    //~ const std::array<double, dim>& v2,
    //~ const std::array<double, dim>& v3)
//~ {
    //~ // Exactamente igual a la original:
    //~ auto lambdas = barycentric_coordinates(p, p0, p1, p2, p3);

    //~ std::array<double, dim> res{};
    //~ for (int d = 0; d < dim; ++d) {
        //~ res[d] = lambdas[0] * v0[d]
               //~ + lambdas[1] * v1[d]
               //~ + lambdas[2] * v2[d]
               //~ + lambdas[3] * v3[d];
    //~ }

    //~ return res;
//~ }

void setVec (double *vec, int size, int dim, double val = 0.0){

for (int v=0;v<size;v++)
    vec[v] = val;
  
}

// void ReMesher::MapNodal(double *vfield, double *o_field){
  // #ifndef REMESH_OMEGA_H
  // MapNodalVectorRaw<3>        (vfield, o_field);
  // #else              
  // MapNodalVector<3> (m_mesh, vfield, o_field); 
  // #endif
// }


void ReMesher::MapNodal(double *vfield, double *o_field){
  const int dim = m_dom->m_dim;
#ifndef REMESH_OMEGA_H
  if(dim==2)
    MapNodalVectorRaw<2>(vfield, o_field);
  else
    MapNodalVectorRaw<3>(vfield, o_field);
#else
  if(dim==2)
  MapNodalVector<2>(m_mesh, vfield, o_field);
  else
  MapNodalVector<3>(m_mesh, vfield, o_field);    
#endif
}


void ReMesher::MapNodalScalar(double *vfield, double *o_field){
  #ifndef REMESH_OMEGA_H
  //MapNodalVectorRaw<1>        (vfield, o_field);
  MapNodalScalarRaw(vfield, o_field);
  #else              

  #endif
  
}

void ReMesher::MapElem(double *vfield, double *o_field, int field_dim){
  #ifndef REMESH_OMEGA_H ////MMG
  MapElemVectorRaw        (vfield, o_field,field_dim);
  //MapElemVectorL2(vfield, o_field,field_dim);
  const int dim = m_dom->m_dim;
  
  #else     
  if (dim ==2)
    MapElemVector<2> (m_mesh, vfield, o_field,field_dim); 
  else
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
  
  bool m_map_momentum = false;
  
  cout << "WRITING DOMAIN "<<m_node_count<<" NODES "<<m_elem_count<<"ELEMS"<<endl;  
  #ifdef REMESH_OMEGA_H
    m_node_count = m_mesh.nverts();
    m_elem_count = m_mesh.nelems();
  #endif

  //memcpy_t(m_->m_elnod, elnod_h, sizeof(int) * dom->m_elem_count * m_dom->m_nodxelem); 
  int dim = m_dom->m_dim;
  
  double *ufield  = new double [dim*m_node_count];   
  //std::vector<double> ufield(3*m_node_count, 0.0);  
  double *vfield   = new double [dim*m_node_count]; 
  double *afield   = new double [dim*m_node_count]; 
  double *pafield   = new double [dim*m_node_count];  //prev_a
  
  double *cforce   = new double [dim*m_node_count];   
  
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
  double *idetF   = new double [m_elem_count];  //Inverse Deformation gradient
  
  double *Tfield  = new double [m_node_count];
  
  double rho_0 = m_dom->rho_0[0];

  //IF MAP MOPMENTUM
  double3 total_momentum_old = make_double3(0.0, 0.0, 0.0);
  for (int i = 0; i < m_dom->m_node_count; i++) {
      total_momentum_old.x += m_dom->m_mdiag[i] * m_dom->v[m_dom->m_dim * i];
      total_momentum_old.y += m_dom->m_mdiag[i] * m_dom->v[m_dom->m_dim * i + 1];
      
      if (dim == 3)
      total_momentum_old.z += m_dom->m_mdiag[i] * m_dom->v[m_dom->m_dim * i + 2];
  }


  // //NEW CONSERVATIVE MOMENTUM MAP
  // double* p_elem_old = new double[m_dom->m_elem_count * 3];
  // double3 pelemold_sum= make_double3(0, 0, 0);
  // for (int e = 0; e < m_dom->m_elem_count; ++e) {
      // double3 p_sum = make_double3(0, 0, 0);
      // for (int a = 0; a < m_dom->m_nodxelem; ++a) {
          // int n = m_dom->m_elnod[e * m_dom->m_nodxelem + a];
          // p_sum.x += m_dom->m_mdiag[n] * m_dom->v[3*n];
          // p_sum.y += m_dom->m_mdiag[n] * m_dom->v[3*n + 1];
          // p_sum.z += m_dom->m_mdiag[n] * m_dom->v[3*n + 2];
      // }
      // p_elem_old[3*e] = p_sum.x / m_dom->m_nodxelem;
      // p_elem_old[3*e + 1] = p_sum.y / m_dom->m_nodxelem;
      // p_elem_old[3*e + 2] = p_sum.z / m_dom->m_nodxelem;
      
      // pelemold_sum.x +=p_elem_old[3*e];
      // pelemold_sum.y+=p_elem_old[3*e+1];
      // pelemold_sum.z+=p_elem_old[3*e+2];
  // }

  
  for (int i=0;i<dim*m_node_count;i++){
    afield[i]=0.0;
    pafield[i]=0.0;
    vfield[i]=0.0;
    cforce[i] = 0.0;
  }

  cout <<"MAP NODAL"<<endl;

  MapNodal(ufield,   m_dom->u); //new , old

  /////// IF MAP VEL DIRECTLY
  
  if (!m_map_momentum){
    MapNodal(vfield,   m_dom->v); //DOES NOT CONS MOMENTUM
  }
  if (m_dom->m_remesh_map_acc){
      MapNodal(afield,   m_dom->a);
      MapNodal(pafield,  m_dom->prev_a);
  }
  
  if (m_dom->m_thermal){
    MapNodalScalar(Tfield,   m_dom->T);
  }

  //MapNodal(pafield,  m_dom->prev_a);
      
  cout <<"DONE. Setting Volumes for "<<m_elem_count<<" elements"<<endl;
  double *volumes=new double[m_elem_count];
  double vol = 0.0;
  for (int i=0;i<m_elem_count;i++){
      int n0 = m_elnod[4 * i];
      int n1 = m_elnod[4 * i + 1];
      int n2 = m_elnod[4 * i + 2];
      int n3 = m_elnod[4 * i + 3];
      int dim = m_dom->m_dim;
      
      if (dim == 3 ){
        double p0 []= {m_x[3*n0], m_x[3*n0+1], m_x[3*n0+2]};
        double p1 [] = {m_x[3*n1], m_x[3*n1+1], m_x[3*n1+2]};
        double p2 []= {m_x[3*n2], m_x[3*n2+1], m_x[3*n2+2]};
        double p3 [] = {m_x[3*n3], m_x[3*n3+1], m_x[3*n3+2]};
        volumes[i]=tet_volume(p0,p1,p2,p3);
      } else { /////// AREA WEIGHT
        double p0[] = {m_x[2*n0], m_x[2*n0+1]};
        double p1[] = {m_x[2*n1], m_x[2*n1+1]};
        double p2[] = {m_x[2*n2], m_x[2*n2+1]};
        double p3[] = {m_x[2*n3], m_x[2*n3+1]};
      // Chequeo de índices de nodos
      if (n0 < 0 || n0 >= m_node_count ||
          n1 < 0 || n1 >= m_node_count ||
          n2 < 0 || n2 >= m_node_count ||
          n3 < 0 || n3 >= m_node_count) {
          std::cerr << "WARNING: Element " << i << " has invalid node indices: "
                    << n0 << ", " << n1 << ", " << n2 << ", " << n3
                    << " (m_node_count=" << m_node_count << ")\n";
      }
        volumes[i] = quad_area(p0, p1, p2, p3);        
        
      }
      vol+=volumes[i];
  }
  cout << "New Element Volume: "<<vol<<endl;
  
  //FindMapElemClosest();
  if (dim == 2)
    FindMapElemClosest<2>();
  else if (dim == 3)
    FindMapElemClosest<3>();
  cout << "foundelem closrsrt"<<endl;
  MapElem(esfield,   m_dom->pl_strain);
  cout << "map pressure"<<endl;
  MapElem(pfield,    m_dom->p);
  
  ///// IF MAP MOMENTUM FOM ELEMENT
  cout << "Map eleme "<<endl;
    double* p_elem_new = new double[m_elem_count * dim];
    //MapElem(p_elem_new, p_elem_old, 3);  // 3 componentes



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
  cout <<"Deleting BCs"<<endl;
  if (m_dom->bc_count[0]>0){free_t(m_dom->bcx_nod);  free_t(m_dom->bcx_val); }
  if (m_dom->bc_count[1]>0){free_t(m_dom->bcy_nod);  free_t(m_dom->bcy_val);  }
  if (m_dom->bc_count[2]>0){free_t(m_dom->bcz_nod);   free_t(m_dom->bcz_val);}
  //~ free_t(bcx_nod);  free_t(bcy_nod);  free_t(bcz_nod);
  //~ free_t(bcx_val);  free_t(bcy_val);  free_t(bcz_val);      
  cout << "Done"<<endl;
  
  m_dom->bc_count[0]=m_dom->bc_count[1]=m_dom->bc_count[2]=0;
  m_dom->bcx_nod_h.clear();  m_dom->bcy_nod_h.clear();  m_dom->bcz_nod_h.clear();
  m_dom->bcx_val_h.clear();  m_dom->bcy_val_h.clear();  m_dom->bcz_val_h.clear();

  
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
    

  memcpy_t(m_dom->u,        ufield, sizeof(double) * m_dom->m_node_count * dim);  
  
  /// THIS SHOULD BE MAINTAINED!!!
  memcpy_t(m_dom->v,        vfield, sizeof(double) * m_dom->m_node_count * dim);
  memcpy_t(m_dom->m_vprev,  vfield, sizeof(double) * m_dom->m_node_count * dim);  
  
  memcpy_t(m_dom->a,        afield, sizeof(double) * m_dom->m_node_count * dim);   
  memcpy_t(m_dom->prev_a,  pafield, sizeof(double) * m_dom->m_node_count * dim);   
  memcpy_t(m_dom->contforce,cforce, sizeof(double) * m_dom->m_node_count * dim);

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
    //THIS IS CRITIC 
    for (int i = 0; i < m_dom->m_elem_count*m_dom->m_nodxelem; ++i)  m_dom->m_dTedt[i] = 0.0; //THERMAL EXPANSION IS BEFORE dTde calac
  }
  
  for (int i=0;i<m_dom->m_node_count*m_dom->m_dim;i++){
    m_dom->m_fe[i]= m_dom->m_fi[i] = 0.0;
  }
   
  m_dom->AssignMatAddress();
  const Material_ *matt  = &m_dom->materials[0];
  cout << "G "<<matt->Elastic().G()<<endl;
  
  ////// AFTER ALL COPIES //////
  cout << "copying"<<endl;
  memcpy_t(m_dom->x,      m_x, m_dom->m_dim*sizeof(double) * m_dom->m_node_count);       
  memcpy_t(m_dom->m_elnod,  m_elnod, 4*sizeof(int) * m_dom->m_elem_count);  
    cout << "Setting conn"<<endl;
  m_dom->setNodElem(m_elnod);     

  //Now calculate new mass, after mapping volume
  cout << "Calculating Volumes"<<endl;
  m_dom->CalcNodalVol();
  cout << "Total Volume "<<m_dom->m_vol_tot<<endl;
  cout << "Calculating mass from Vols"<<endl;
  m_dom->CalcNodalMassFromVol();
  for (int i=0;i<m_node_count;i++){ //Or already dom_d->m_node_count since domain changed
    if (m_dom->m_mdiag[i]<1.0e-10)
      cout << "ERROR, SMALL MASS ON NODE "<<i<<endl;
  }
  //if (m_dom->m_remesh_map_vel)
 cout << "recovering velocities"<<endl;
  if (m_map_momentum){
    for (int i=0;i<m_node_count;i++){ //Or already dom_d->m_node_count since domain changed
      for (int d=0;d<m_dom->m_dim;d++)
        vfield[m_dom->m_dim*i+d] /=m_dom->m_mdiag[i];
    }


    /// --- 3. Calcular momento DESPUÉS del mapeo ---
    double3 total_momentum_new = make_double3(0.0, 0.0, 0.0);
    for (int i = 0; i < m_node_count; i++) {
        total_momentum_new.x += m_dom->m_mdiag[i] * vfield[m_dom->m_dim * i];
        total_momentum_new.y += m_dom->m_mdiag[i] * vfield[m_dom->m_dim * i + 1];
        if (dim == 3)total_momentum_new.z += m_dom->m_mdiag[i] * vfield[m_dom->m_dim * i + 2];
    }

    // --- 4. Aplicar corrección si hay discrepancia (>1% error) ---
    double3 correction_factor = make_double3(1.0, 1.0, 1.0);
    if (fabs(total_momentum_old.x) > 1e-10) correction_factor.x = total_momentum_old.x / total_momentum_new.x;
    if (fabs(total_momentum_old.y) > 1e-10) correction_factor.y = total_momentum_old.y / total_momentum_new.y;
    if (dim == 3)
      if (fabs(total_momentum_old.z) > 1e-10) correction_factor.z = total_momentum_old.z / total_momentum_new.z;

    
    // --- 5. Verificación final (opcional, para debug) ---
    double3 final_momentum = make_double3(0.0, 0.0, 0.0);
    for (int i = 0; i < m_node_count; i++) {
        final_momentum.x += m_dom->m_mdiag[i] * vfield[m_dom->m_dim * i];
        final_momentum.y += m_dom->m_mdiag[i] * vfield[m_dom->m_dim * i + 1];
        if (dim == 3)
        final_momentum.z += m_dom->m_mdiag[i] * vfield[m_dom->m_dim * i + 2];
    }

    memcpy_t(m_dom->m_vprev,  vfield, sizeof(double) * m_dom->m_node_count * dim); 
    memcpy_t(m_dom->v,        vfield, sizeof(double) * m_dom->m_node_count * dim); 
        
    cout << "Momentum:\n";
    cout << " - Before:  (" << total_momentum_old.x << ", " << total_momentum_old.y << ", " << total_momentum_old.z << ")\n";
    cout << " - After: (" << final_momentum.x << ", " << final_momentum.y << ", " << final_momentum.z << ")\n";
      
  
  }/////// IF MAP MOMENTUM
  
  //// RECALCULATED FROM MOMENTUM
  ///////if (m_dom->m_remesh_map_vel)
  /////memcpy_t(m_dom->v,        vfield, sizeof(double) * m_dom->m_node_count * 3); 
          
  
  
  //~ ReMapBCs(bcx_nod,bcx_val,m_dom->bcx_nod, m_dom->bcx_val, bccount[0]);
  //~ ReMapBCs(bcy_nod,bcy_val,m_dom->bcy_nod, m_dom->bcy_val, bccount[1]);
  //~ ReMapBCs(bcz_nod,bcz_val,m_dom->bcz_nod, m_dom->bcz_val, bccount[2]);

  
  cout << "deleting "<<endl;
  delete[] ufield; 
  delete[] vfield;
  delete[] afield;
  delete[] pafield;
  delete[] esfield;
  delete[] pfield;
  delete[] sigfield;
  delete[] syfield;
  delete[] psfield;
  delete[] str_rate;
  delete[] rot_rate;
  delete[] tau;
  delete[] rho;
  delete[] vol_0;
  delete[] idetF;
  //delete [] bcx_nod,bcy_nod,bcz_nod,bcx_val,bcy_val,bcz_val;
  delete [] Tfield;
  cout << "MESH CHANGED"<<endl;


  //delete[] p_elem_old;
  delete[] p_elem_new;
  delete[] volumes;
  //AFTER MAP
  //THIS CRASHES
  //free_t(m_closest_elem);
  

}


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

// int ReMesher::find_closest_node(const std::array<double, 3>& x) {
    // int closest = -1;
    // double min_dist_sq = 1.0e20;

    // for (int i = 0; i < m_dom->m_node_count; ++i) {
        // double dx = x[0] - m_dom->x[3 * i];
        // double dy = x[1] - m_dom->x[3 * i + 1];
        // double dz = x[2] - m_dom->x[3 * i + 2];
        // double dist_sq = dx * dx + dy * dy + dz * dz;

        // if (dist_sq < min_dist_sq) {
            // min_dist_sq = dist_sq;
            // closest = i;
        // }
    // }
    // return closest;
// }

template<int dim>
int ReMesher::find_closest_node(const std::array<double, dim>& x)
{
    static_assert(dim == 2 || dim == 3, "dim must be 2 or 3");

    int closest = -1;
    double min_dist_sq = 1.0e20;

    for (int i = 0; i < m_dom->m_node_count; ++i) {
        double dist_sq = 0.0;

        #pragma unroll
        for (int d = 0; d < dim; ++d) {
            double dx = x[d] - m_dom->x[dim * i + d];
            dist_sq += dx * dx;
        }

        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            closest = i;
        }
    }

    return closest;
}

///////////////NEW, BOTH 2D and 3D
template <int dim>
void ReMesher::MapNodalVectorRaw(double *vfield, double *o_field)
{
    int notfound = 0;
    int notsame  = 0;

    cout << "MAP NODAL VECTOR RAW (MMG), dim=" << dim << endl;

    for (int vert = 0; vert < m_node_count; vert++) {

        // Inicializo el vector de destino
        for (int d = 0; d < dim; d++)
            vfield[dim * vert + d] = 0.0;

        // Coordenadas del nodo actual
        std::array<double, dim> x;
        for (int d = 0; d < dim; d++)
            x[d] = m_x[dim * vert + d];

        bool found_samenode = false;
        const double tol = 1e-6;

        // ---- SAME NODE ----
        for (int v = 0; v < m_dom->m_node_count; v++) {

            double dist = 0.0;
            for (int d = 0; d < dim; d++) {
                double dx = x[d] - m_dom->x[dim * v + d];
                dist += dx * dx;
            }

            if (dist < tol) {
                for (int d = 0; d < dim; d++)
                    vfield[dim * vert + d] = o_field[dim * v + d];
                found_samenode = true;
                break;
            }
        }

        if (found_samenode) continue;

        notsame++;
        bool found = false;

        // ---- ELEMENT SEARCH ----
        for (int e = 0; e < m_dom->m_elem_count; e++) {

            int n[4];
            for (int i = 0; i < 4; i++)
                n[i] = m_dom->m_elnod[4 * e + i];

            // Creo puntos en 3D siempre, z=0 si es 2D
            std::array<double, 3> p3[4];
            for (int i = 0; i < 4; i++) {
                p3[i][0] = m_dom->x[dim * n[i] + 0];
                p3[i][1] = m_dom->x[dim * n[i] + 1];
                p3[i][2] = (dim == 2 ? 0.0 : m_dom->x[dim * n[i] + 2]);
            }

            std::array<double, 3> x3 = {x[0], x[1], (dim == 2 ? 0.0 : x[2])};

            // Calcula lambdas en 3D
            std::array<double, 4> w = stable_barycentric(x3, p3[0], p3[1], p3[2], p3[3]);

            // Interpolación vectorial
            std::array<double, dim> interp{};
            for (int i = 0; i < 4; i++)
                for (int d = 0; d < dim; d++)
                    interp[d] += w[i] * o_field[dim * n[i] + d];

            for (int d = 0; d < dim; d++)
                vfield[dim * vert + d] = interp[d];

            found = true;
            break;
        }

        // ---- FALLBACK ----
        if (!found) {
            int closest = find_closest_node<dim>(x);
            for (int d = 0; d < dim; d++)
                vfield[dim * vert + d] = o_field[dim * closest + d];
            notfound++;
        }
    }

    cout << "Not inside element: " << notfound
         << " (" << 100.0 * notfound / m_node_count << "%)\n";
    cout << "Not same node: " << notsame
         << " (" << 100.0 * notsame / m_node_count << "%)\n";
}

void ReMesher::MapNodalScalarRaw(double* sfield, double* o_field)
{
    cout << "MAPPING NODAL SCALAR"<<endl;
    for (int vert = 0; vert < m_node_count; vert++) {

        double x[3];
        for (int d = 0; d < 3; d++)
            x[d] = m_x[3*vert + d];

        bool found = false;

        for (int i = 0; i < m_dom->m_elem_count; i++) {

            int n0 = m_dom->m_elnod[4*i+0];
            int n1 = m_dom->m_elnod[4*i+1];
            int n2 = m_dom->m_elnod[4*i+2];
            int n3 = m_dom->m_elnod[4*i+3];

            std::array<double, 3> p0 = {m_dom->x[3*n0], m_dom->x[3*n0+1], m_dom->x[3*n0+2]};
            std::array<double, 3> p1 = {m_dom->x[3*n1], m_dom->x[3*n1+1], m_dom->x[3*n1+2]};
            std::array<double, 3> p2 = {m_dom->x[3*n2], m_dom->x[3*n2+1], m_dom->x[3*n2+2]};
            std::array<double, 3> p3 = {m_dom->x[3*n3], m_dom->x[3*n3+1], m_dom->x[3*n3+2]};

            auto lambda = stable_barycentric({x[0],x[1],x[2]}, p0,p1,p2,p3);
            const double eps = 1e-10;

            double sum =
                lambda[0] + lambda[1] + lambda[2] + lambda[3];

            bool inside =
                lambda[0] >= -eps && lambda[1] >= -eps &&
                lambda[2] >= -eps && lambda[3] >= -eps &&
                fabs(sum - 1.0) < 1e-8; //Critic, check for garbage
            if (inside) {

                sfield[vert] =
                      lambda[0]*o_field[n0]
                    + lambda[1]*o_field[n1]
                    + lambda[2]*o_field[n2]
                    + lambda[3]*o_field[n3];

                double omin = std::min({o_field[n0], o_field[n1],
                                        o_field[n2], o_field[n3]});
                double omax = std::max({o_field[n0], o_field[n1],
                                        o_field[n2], o_field[n3]});

                if (sfield[vert] < omin - 1e-8 || sfield[vert] > omax + 1e-8) {
                    cout << "WARNING: non-monotone interp at vert "
                         << vert << endl;
                }

                found = true;
                break;
            }
        }

        if (!found) {
            int nn = find_closest_node<3>({x[0],x[1],x[2]});

            if (nn < 0 || !std::isfinite(o_field[nn])) {
                sfield[vert] = 0.0;
                cout << "WARNING: fallback failed at vert "
                     << vert << endl;
            } else {
                sfield[vert] = o_field[nn];
            }
            
            //sfield[vert] = o_field[nn];

        }
    }
    
    double vmax = 0.0;
    double vmin = 1.0e6;
    for (int vert=0;vert<m_node_count;vert++){
      if (sfield[vert]>vmax) vmax = sfield[vert];
      if (sfield[vert]<vmin) vmin = sfield[vert];
    }
    cout <<"MAX INTERPOLATED VALUE "<< vmax<<"; MIN "<<vmin<<endl;
    
     vmax = 0.0;
     vmin = 1.0e6;
    for (int vert=0;vert<m_dom->m_node_count;vert++){
      if (o_field[vert]>vmax) vmax = o_field[vert];
      if (o_field[vert]<vmin) vmin = o_field[vert];
    }
    cout <<"MAX ORIGINAL VALUE "<< vmax<<"; MIN "<<vmin<<endl;    

}

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



//~ void ReMesher::FindMapElemClosest() {
  //~ int notinside = 0;
    //~ for (int elem = 0; elem < m_elem_count; elem++) {
        //~ std::array<double, 3> barycenter = {0.0, 0.0, 0.0};
        
        //~ //////Calculate barycenter of current new element
        //~ for (int en = 0; en < 4; en++) {
            //~ int v = m_elnod[elem * 4 + en];
            //~ barycenter[0] += m_x[3 * v];
            //~ barycenter[1] += m_x[3 * v + 1];
            //~ barycenter[2] += m_x[3 * v + 2];
        //~ }
        //~ barycenter[0] /= 4.0;
        //~ barycenter[1] /= 4.0;
        //~ barycenter[2] /= 4.0;

        //~ bool found_containing = false;
        //~ int closest_elem = -1;
        //~ double min_distance = std::numeric_limits<double>::max();
        //~ std::array<double, 3> closest_barycenter;

        //~ ////////First pass: Try to find an element that contains the barycenter
        //~ for (int i = 0; i < m_dom->m_elem_count; i++) {
            //~ int n0 = m_dom->m_elnod[4 * i];
            //~ int n1 = m_dom->m_elnod[4 * i + 1];
            //~ int n2 = m_dom->m_elnod[4 * i + 2];
            //~ int n3 = m_dom->m_elnod[4 * i + 3];

            //~ std::array<double, 3> p0 = {m_dom->x[3 * n0], m_dom->x[3 * n0 + 1], m_dom->x[3 * n0 + 2]};
            //~ std::array<double, 3> p1 = {m_dom->x[3 * n1], m_dom->x[3 * n1 + 1], m_dom->x[3 * n1 + 2]};
            //~ std::array<double, 3> p2 = {m_dom->x[3 * n2], m_dom->x[3 * n2 + 1], m_dom->x[3 * n2 + 2]};
            //~ std::array<double, 3> p3 = {m_dom->x[3 * n3], m_dom->x[3 * n3 + 1], m_dom->x[3 * n3 + 2]};

            //~ ////////Reuse your existing barycentric coordinate function
            //~ std::array<double, 4> lambdas = stable_barycentric(barycenter, p0, p1, p2, p3);

            //~ if (lambdas[0] >= -1.0e-8 && lambdas[1] >= -1.0e-8 && 
                //~ lambdas[2] >= -1.0e-8 && lambdas[3] >= -1.0e-8) {
                //~ m_closest_elem[elem] = i;
                //~ found_containing = true;
                
                //~ ////////Diagnostic output
                //~ if (elem < 5) { ////////Print first few elements for verification
                    //~ std::cout << "Element " << elem << " barycenter INSIDE old element " << i 
                              //~ << " with barycentric coords: (" 
                              //~ << lambdas[0] << ", " << lambdas[1] << ", " 
                              //~ << lambdas[2] << ", " << lambdas[3] << ")\n";
                //~ }
                //~ break;
            //~ }

            //~ /////////////////If not inside, compute distance for fallback
            //~ std::array<double, 3> old_barycenter = {
                //~ (p0[0] + p1[0] + p2[0] + p3[0]) / 4.0,
                //~ (p0[1] + p1[1] + p2[1] + p3[1]) / 4.0,
                //~ (p0[2] + p1[2] + p2[2] + p3[2]) / 4.0
            //~ };

            //~ double distance = std::pow(barycenter[0] - old_barycenter[0], 2) +
                             //~ std::pow(barycenter[1] - old_barycenter[1], 2) +
                             //~ std::pow(barycenter[2] - old_barycenter[2], 2);

            //~ if (distance < min_distance) {
                //~ min_distance = distance;
                //~ closest_elem = i;
                //~ closest_barycenter = old_barycenter;
            //~ }
        //~ }

        //~ ///////////////Fallback to closest element if no containing element found
        //~ if (!found_containing) {
            //~ m_closest_elem[elem] = closest_elem;
            
            //~ ////////////////Diagnostic output
            //~ if (elem < 5) {
                //~ std::cout << "Element " << elem << " barycenter NOT INSIDE any element. "
                          //~ << "Using closest element " << closest_elem 
                          //~ << " with distance " << std::sqrt(min_distance) << "\n";
                //~ std::cout << "New barycenter: (" << barycenter[0] << ", " 
                          //~ << barycenter[1] << ", " << barycenter[2] << ")\n";
                //~ std::cout << "Old barycenter: (" << closest_barycenter[0] << ", " 
                          //~ << closest_barycenter[1] << ", " << closest_barycenter[2] << ")\n";
            //~ }
        //~ } else {
          
          //~ notinside++;
        //~ }
    //~ }//ELEM
    //~ cout << "Not Inside  Elements : " << notinside<<"( "<< notinside/m_elem_count*100.0<<"%)"<<endl;

//~ }

template <int dim>
void ReMesher::FindMapElemClosest() {
    int notinside = 0;

    constexpr int nodes_per_elem = 4; // siempre 4 nodos (quads en 2D, tetras en 3D)

    for (int elem = 0; elem < m_elem_count; elem++) {

        /////////// Calcular barycenter del elemento actual
        std::array<double, 3> barycenter = {0.0, 0.0, 0.0};
        for (int en = 0; en < nodes_per_elem; en++) {
            int v = m_elnod[elem * nodes_per_elem + en];
            barycenter[0] += m_x[dim * v + 0];
            barycenter[1] += m_x[dim * v + 1];
            barycenter[2] += (dim == 3) ? m_x[dim * v + 2] : 0.0;
        }
        for (int d = 0; d < 3; d++)
            barycenter[d] /= static_cast<double>(nodes_per_elem);

        bool found_containing = false;
        int closest_elem = -1;
        double min_distance = std::numeric_limits<double>::max();
        std::array<double, 3> closest_barycenter{};

        /////////// Primer pase: buscar elemento que contenga el barycenter
        for (int i = 0; i < m_dom->m_elem_count; i++) {

            // Preparar nodos para stable_barycentric
            std::array<double, 3> p0, p1, p2, p3;
            int n0 = m_dom->m_elnod[4*i + 0];
            int n1 = m_dom->m_elnod[4*i + 1];
            int n2 = m_dom->m_elnod[4*i + 2];
            int n3 = m_dom->m_elnod[4*i + 3];

            p0[0] = m_dom->x[dim*n0 + 0]; p0[1] = m_dom->x[dim*n0 + 1]; p0[2] = (dim==3) ? m_dom->x[dim*n0 + 2] : 0.0;
            p1[0] = m_dom->x[dim*n1 + 0]; p1[1] = m_dom->x[dim*n1 + 1]; p1[2] = (dim==3) ? m_dom->x[dim*n1 + 2] : 0.0;
            p2[0] = m_dom->x[dim*n2 + 0]; p2[1] = m_dom->x[dim*n2 + 1]; p2[2] = (dim==3) ? m_dom->x[dim*n2 + 2] : 0.0;
            p3[0] = m_dom->x[dim*n3 + 0]; p3[1] = m_dom->x[dim*n3 + 1]; p3[2] = (dim==3) ? m_dom->x[dim*n3 + 2] : 0.0;

            // Llamada a la función original
            std::array<double,4> lambdas = stable_barycentric(barycenter, p0, p1, p2, p3);

            bool inside = true;
            for (int n = 0; n < nodes_per_elem; n++)
                if (lambdas[n] < -1.0e-8) inside = false;

            if (inside) {
                m_closest_elem[elem] = i;
                found_containing = true;

                if (elem < 5) {
                    std::cout << "Element " << elem << " barycenter INSIDE old element " << i 
                              << " with barycentric coords: ";
                    for (int n = 0; n < nodes_per_elem; n++)
                        std::cout << lambdas[n] << " ";
                    std::cout << "\n";
                }
                break;
            }

            ////////// Si no está dentro, calcular distancia para fallback
            std::array<double, 3> old_barycenter = {0.0, 0.0, 0.0};
            for (int d = 0; d < 3; d++)
                old_barycenter[d] = (p0[d] + p1[d] + p2[d] + p3[d]) / 4.0;

            double distance = 0.0;
            for (int d = 0; d < dim; d++)
                distance += (barycenter[d] - old_barycenter[d]) * (barycenter[d] - old_barycenter[d]);

            if (distance < min_distance) {
                min_distance = distance;
                closest_elem = i;
                closest_barycenter = old_barycenter;
            }
        }

        /////////// Fallback si no encuentra elemento contenedor
        if (!found_containing) {
            m_closest_elem[elem] = closest_elem;

            if (elem < 5) {
                std::cout << "Element " << elem << " barycenter NOT INSIDE any element. "
                          << "Using closest element " << closest_elem
                          << " with distance " << std::sqrt(min_distance) << "\n";
                std::cout << "New barycenter: ";
                for (int d = 0; d < dim; d++) std::cout << barycenter[d] << " ";
                std::cout << "\nOld barycenter: ";
                for (int d = 0; d < dim; d++) std::cout << closest_barycenter[d] << " ";
                std::cout << "\n";
            }
        } else {
            notinside++;
        }

    } // elem

    std::cout << "Not Inside Elements : " << notinside
              << " ( " << notinside * 100.0 / m_elem_count << "%)\n";
}


void ReMesher::ReMapBCs(int  *old_bc_nod,
                      double *old_bc_val,

                    int  *new_bc_nod,
                    double *new_bc_val,
                    int bc_count) {

}




void ReMesher::ReMapBCsByFace(int* old_bc_nod,
                        double* old_bc_val,
                        int* new_bc_nod,
                        double* new_bc_val,
                        int bc_count) {
  
  
}


//~ //////////////////////////// NEW MIN SQUARES ELEMENT MAP
//~ #include <Eigen/Dense>
//~ #include <Eigen/QR> // Para ColPivHouseholderQR
//~ #include <vector>
//~ #include <array>
//~ #include <cmath>

//~ void ReMesher::MapElemVectorL2(double* new_field, double* old_field, int field_dim) {
    //~ // 1. Precomputar centroides de la malla OLD (si no se tiene ya)
    //~ std::vector<std::array<double, 3>> old_centroids(m_dom->m_elem_count);
    //~ for (int old_elem = 0; old_elem < m_dom->m_elem_count; old_elem++) {
        //~ old_centroids[old_elem] = get_old_elem_centroid(old_elem);
    //~ }

    //~ // 2. Para cada elemento en la NUEVA malla
    //~ for (int new_elem = 0; new_elem < m_elem_count; new_elem++) {
        
        //~ // 3. Obtener centroide del elemento nuevo
        //~ std::array<double, 3> new_centroid = get_new_elem_centroid(new_elem);
        
        //~ // 4. Encontrar vecinos más cercanos (implementación temporal - optimizar con KD-tree después)
        //~ int num_neighbors = std::min(20, m_dom->m_elem_count);
        //~ std::vector<int> neighbor_old_elems = find_nearest_old_elems(new_centroid, num_neighbors, old_centroids);
        
        //~ int n = neighbor_old_elems.size();
        
        //~ // 5. Para cada componente del campo vectorial
        //~ for (int comp = 0; comp < field_dim; comp++) {
            //~ double new_value = 0.0;

            //~ if (n < 4) { // Fallback: promedio simple
                //~ double sum = 0.0;
                //~ for (int old_elem_id : neighbor_old_elems) {
                    //~ sum += old_field[old_elem_id * field_dim + comp];
                //~ }
                //~ new_value = (n > 0) ? sum / n : 0.0;
            //~ } else {
                //~ // 6. Configurar problema de mínimos cuadrados A * x = b
                //~ Eigen::MatrixXd A(n, 4);
                //~ Eigen::VectorXd b(n);

                //~ for (int i = 0; i < n; i++) {
                    //~ int old_elem_id = neighbor_old_elems[i];
                    //~ std::array<double, 3>& c = old_centroids[old_elem_id];
                    
                    //~ A(i, 0) = 1.0;
                    //~ A(i, 1) = c[0];
                    //~ A(i, 2) = c[1];
                    //~ A(i, 3) = c[2];
                    
                    //~ b(i) = old_field[old_elem_id * field_dim + comp];
                //~ }

                //~ // 7. Resolver con QR pivoteado (robusto)
                //~ Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
                //~ if (qr.info() != Eigen::Success) {
                    //~ // Fallback si la descomposición falla
                    //~ double sum = 0.0;
                    //~ for (int old_elem_id : neighbor_old_elems) {
                        //~ sum += old_field[old_elem_id * field_dim + comp];
                    //~ }
                    //~ new_value = sum / n;
                    //~ continue;
                //~ }
                
                //~ Eigen::Vector4d x = qr.solve(b);
                //~ new_value = x[0] + x[1]*new_centroid[0] + x[2]*new_centroid[1] + x[3]*new_centroid[2];
            //~ }
            
            //~ // 8. Asignar valor proyectado
            //~ new_field[new_elem * field_dim + comp] = new_value;
        //~ }
    //~ }
//~ }

// Función auxiliar para encontrar elementos cercanos (versión temporal)
std::vector<int> ReMesher::find_nearest_old_elems(const std::array<double, 3>& point, 
                                                 int num_neighbors,
                                                 const std::vector<std::array<double, 3>>& old_centroids) {
    std::vector<std::pair<double, int>> distances;
    
    for (int i = 0; i < old_centroids.size(); i++) {
        double dist_sq = 0.0;
        for (int d = 0; d < 3; d++) {
            dist_sq += std::pow(point[d] - old_centroids[i][d], 2);
        }
        distances.emplace_back(dist_sq, i);
    }
    
    // Ordenar por distancia
    std::sort(distances.begin(), distances.end(), 
             [](const auto& a, const auto& b) { return a.first < b.first; });
    
    // Tomar los n más cercanos
    std::vector<int> neighbors;
    int n = std::min(num_neighbors, (int)distances.size());
    for (int i = 0; i < n; i++) {
        neighbors.push_back(distances[i].second);
    }
    
    return neighbors;
}

// Función para obtener centroide de elemento OLD
std::array<double, 3> ReMesher::get_old_elem_centroid(int old_elem) {
    std::array<double, 3> centroid = {0.0, 0.0, 0.0};
    for (int en = 0; en < 4; en++) {
        int node_id = m_dom->m_elnod[old_elem * 4 + en];
        centroid[0] += m_dom->x[3 * node_id];
        centroid[1] += m_dom->x[3 * node_id + 1];
        centroid[2] += m_dom->x[3 * node_id + 2];
    }
    centroid[0] /= 4.0;
    centroid[1] /= 4.0;
    centroid[2] /= 4.0;
    return centroid;
}

// Función para obtener centroide de elemento NEW
std::array<double, 3> ReMesher::get_new_elem_centroid(int new_elem) {
    std::array<double, 3> centroid = {0.0, 0.0, 0.0};
    for (int en = 0; en < 4; en++) {
        int node_id = m_elnod[new_elem * 4 + en];
        centroid[0] += m_x[3 * node_id];
        centroid[1] += m_x[3 * node_id + 1];
        centroid[2] += m_x[3 * node_id + 2];
    }
    centroid[0] /= 4.0;
    centroid[1] /= 4.0;
    centroid[2] /= 4.0;
    return centroid;
}

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <iostream>

void ReMesher::MapElemVectorL2(double* new_field, double* old_field, int field_dim) {
    // Precomputar centroides OLD
    std::vector<std::array<double, 3>> old_centroids(m_dom->m_elem_count);
    for (int old_elem = 0; old_elem < m_dom->m_elem_count; old_elem++) {
        old_centroids[old_elem] = get_old_elem_centroid(old_elem);
    }

    for (int new_elem = 0; new_elem < m_elem_count; new_elem++) {
        std::array<double, 3> new_centroid = get_new_elem_centroid(new_elem);
        int num_neighbors = std::min(20, m_dom->m_elem_count);
        std::vector<int> neighbor_old_elems = find_nearest_old_elems(new_centroid, num_neighbors, old_centroids);
        
        int n = neighbor_old_elems.size();
        
        for (int comp = 0; comp < field_dim; comp++) {
            if (n < 4) {
                // Fallback: promedio ponderado por distancia inversa
                double sum_val = 0.0;
                double sum_weight = 0.0;
                
                for (int old_elem_id : neighbor_old_elems) {
                    double dist = distance(new_centroid, old_centroids[old_elem_id]);
                    double weight = 1.0 / (dist * dist + 1e-15);
                    sum_val += weight * old_field[old_elem_id * field_dim + comp];
                    sum_weight += weight;
                }
                new_field[new_elem * field_dim + comp] = (sum_weight > 0) ? sum_val / sum_weight : 0.0;
            } else {
                // Implementación manual de mínimos cuadrados
                double ATA[4][4] = {0};
                double ATb[4] = {0};
                
                // Construir A^T A y A^T b
                for (int i = 0; i < n; i++) {
                    int old_elem_id = neighbor_old_elems[i];
                    auto& c = old_centroids[old_elem_id];
                    double val = old_field[old_elem_id * field_dim + comp];
                    
                    double row[4] = {1.0, c[0], c[1], c[2]};
                    
                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 4; k++) {
                            ATA[j][k] += row[j] * row[k];
                        }
                        ATb[j] += row[j] * val;
                    }
                }
                
                // Resolver sistema normal A^T A x = A^T b
                double x[4];
                if (solve_4x4_system(ATA, ATb, x)) {
                    new_field[new_elem * field_dim + comp] = x[0] + x[1]*new_centroid[0] + 
                                                           x[2]*new_centroid[1] + x[3]*new_centroid[2];
                } else {
                    // Fallback a interpolación ponderada
                    double sum_val = 0.0;
                    double sum_weight = 0.0;
                    for (int old_elem_id : neighbor_old_elems) {
                        double dist = distance(new_centroid, old_centroids[old_elem_id]);
                        double weight = 1.0 / (dist * dist + 1e-15);
                        sum_val += weight * old_field[old_elem_id * field_dim + comp];
                        sum_weight += weight;
                    }
                    new_field[new_elem * field_dim + comp] = (sum_weight > 0) ? sum_val / sum_weight : 0.0;
                }
            }
        }
    }
}

bool ReMesher::solve_4x4_system(double A[4][4], double b[4], double x[4]) {
    // Implementación robusta de eliminación gaussiana con pivoteo
    int pivot[4] = {0, 1, 2, 3};
    
    // Eliminación hacia adelante
    for (int col = 0; col < 4; col++) {
        // Búsqueda de pivote
        int max_row = col;
        double max_val = std::abs(A[col][col]);
        for (int row = col + 1; row < 4; row++) {
            if (std::abs(A[row][col]) > max_val) {
                max_val = std::abs(A[row][col]);
                max_row = row;
            }
        }
        
        if (max_val < 1e-12) {
            return false; // Matriz singular
        }
        
        if (max_row != col) {
            // Intercambiar filas
            for (int i = 0; i < 4; i++) {
                std::swap(A[col][i], A[max_row][i]);
            }
            std::swap(b[col], b[max_row]);
            std::swap(pivot[col], pivot[max_row]);
        }
        
        // Eliminación
        for (int row = col + 1; row < 4; row++) {
            double factor = A[row][col] / A[col][col];
            for (int i = col + 1; i < 4; i++) {
                A[row][i] -= factor * A[col][i];
            }
            b[row] -= factor * b[col];
            A[row][col] = 0.0;
        }
    }
    
    // Sustitución hacia atrás
    for (int row = 3; row >= 0; row--) {
        x[row] = b[row];
        for (int col = row + 1; col < 4; col++) {
            x[row] -= A[row][col] * x[col];
        }
        x[row] /= A[row][row];
    }
    
    // Reordenar solución según pivoteo
    double temp[4];
    for (int i = 0; i < 4; i++) temp[i] = x[i];
    for (int i = 0; i < 4; i++) x[pivot[i]] = temp[i];
    
    return true;
}

// Función de distancia euclidiana
double ReMesher::distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    double dist_sq = 0.0;
    for (int i = 0; i < 3; i++) {
        double diff = a[i] - b[i];
        dist_sq += diff * diff;
    }
    return std::sqrt(dist_sq);
}
  
};
