/*************************************************************************/
/*  Domain_d.C                                                   */
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



#include "Domain_d.h"
#include <iostream>
#include <vector>

#include "Matrix.h"
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

#include "Tensor3.C"
#include "../lib/LSDynaReader/src/lsdynaReader.h"

#include "parallel_for_each.h"

#include "Mesh.h"

using namespace std;
using namespace LS_Dyna;

namespace MetFEM {

dev_t bool areFacesEqual(const Face& f1, const Face& f2) {
    if (f1.n_nodes != f2.n_nodes) return false;

    bool matched[MAX_FACE_NODES] = {false};

    for (int i = 0; i < f1.n_nodes; i++) {
        bool found = false;
        for (int j = 0; j < f2.n_nodes; j++) {
            if (!matched[j] && f1.nodes[i] == f2.nodes[j]) {
                matched[j] = true;
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true; // todos los nodos coinciden (orden ignorado)
}

dev_t void addFace(Face faceList[], int& faceCount, const Face& newFace) {
    for (int i = 0; i < faceCount; i++) {
        if (areFacesEqual(faceList[i], newFace)) {
            faceList[i].count++;
            if (faceList[i].count == 2) {
                faceList[i].other_elem = newFace.elem_id;
            }
            return;
        }
    }
    faceList[faceCount] = newFace;
    faceList[faceCount].count = 1;
    faceList[faceCount].other_elem = -1;
    faceCount++;
}

dev_t void addElementFaces(Face faceList[], int& faceCount,
                                const int element[], int elem_id,
                                int ELFAC, int FACENOD,
                                const int local_faces[][MAX_FACE_NODES]) 
{
    for (int f = 0; f < ELFAC; f++) {
        Face newFace;
        newFace.n_nodes = FACENOD;
        newFace.elem_id = elem_id;
        newFace.count = 0;
        newFace.other_elem = -1;

        for (int n = 0; n < FACENOD; n++) {
            newFace.nodes[n] = element[ local_faces[f][n] ]; /////CHANGE TO CURRENT IF MIXED ELEMENT
        }

        addFace(faceList, faceCount, newFace);
    }
}

dev_t void Domain_d::SearchExtNodes() {
    
    if (m_dim == 2)
      set2DFacesValues();
    
    #ifndef CUDA_BUILD
    int elements[ELNOD]; 
    #else 
    int elements[4];
    #endif
    
    m_faceCount = 0;
    
    // Extraer caras/aristas de cada elemento
    for (int i = 0; i < m_elem_count; i++) {
        for (int ne = 0; ne < m_nodxelem; ne++)
            elements[ne] = m_elnod[m_nodxelem * i + ne];

        addElementFaces(faceList, m_faceCount,
                        elements, i,
                        ELFAC, FACENOD,
                        local_faces);
    }
    printf("done. Face count: %d\n", m_faceCount);

    // Array de nodos externos
    for (int n = 0; n < m_node_count; n++)
        ext_nodes[n] = false;
    ext_nodes_count = 0;

    // Inicializar listas de vecinos
    for (int i = 0; i < m_elem_count; i++) {
        for (int j = 0; j < ELFAC; j++)
            m_elem_neigh[ELFAC * i + j] = -1;  // -1 = sin vecino
        m_elem_neigh_count[i] = 0;
    }

    // Identificar nodos externos y vecinos internos
    int ext_faces = 0;
    for (int i = 0; i < m_faceCount; i++) {
        if (faceList[i].count == 1) { // Cara externa
            //printf("External Face %d nodes: ", ext_faces);
            for (int j = 0; j < FACENOD; j++) {
                ext_nodes[faceList[i].nodes[j]] = true;
                //printf("%d ", faceList[i].nodes[j]); // LISTA DE NODOS//DEBUG
            }
            ext_faces++;
            //printf("\n");//DEBUG
        }
        else if (faceList[i].count == 2) { // Cara interna
            int e1 = faceList[i].elem_id;
            int e2 = faceList[i].other_elem;

            // Agregar e2 a lista de vecinos de e1
            m_elem_neigh[ELFAC * e1 + m_elem_neigh_count[e1]] = e2;
            m_elem_neigh_count[e1]++;

            // Agregar e1 a lista de vecinos de e2
            m_elem_neigh[ELFAC * e2 + m_elem_neigh_count[e2]] = e1;
            m_elem_neigh_count[e2]++;
        }
    }

    // Contar nodos externos
    for (int i = 0; i < m_node_count; i++) {
        if (ext_nodes[i]) {
            ext_nodes_count++;
        }
    }

    // Calcular áreas en nodos externos
    CalcExtFaceAreas();

    double area = 0.0;
    for (int i = 0; i < m_node_count; i++) {
        if (ext_nodes[i]) {
            area += node_area[i];
        }
    }
    printf("Total External Nodal Area: %.4e, external nodes: %d\n", area,ext_nodes_count);
}



///////////////////////////////////////////////////

dev_t void Domain_d::CalcExtFaceAreas(){
  
  int top_faces = 0;
  double top_area=0.0;
  #ifndef CUDA_BUILD
  bool elem_flags[m_elem_count];
  #else
  bool *elem_flags;
  malloc_t(elem_flags, bool, m_elem_count);
  #endif
      
  for (int e=0;e<m_elem_count;e++)elem_flags[e]=false;

  for (int i = 0; i < m_node_count; i++)
      node_area[i] = 0.0;
  for (int i = 0; i < m_elem_count; i++)
    m_elem_area[i] = 0.0;
        
  if (m_dim == 2 && FACENOD == 2) {

    for (int i = 0; i < m_faceCount; i++) {

        if (faceList[i].count == 1) {   // External edge (2D)

            int n0 = faceList[i].nodes[0];
            int n1 = faceList[i].nodes[1];

            double2 p0 = getPosVec2(n0);
            double2 p1 = getPosVec2(n1);

            double dx = p1.x - p0.x;
            double dy = p1.y - p0.y;
            double length = sqrt(dx*dx + dy*dy);

            // --- 2D PLANO ---
            double edge_measure = length;

            // --- 2D AXISIMÉTRICO (si aplica) ---
            // double rmid = 0.5 * (p0.x + p1.x);
            // double edge_measure = 2.0 * M_PI * rmid * length;

            // Distribuir a nodos
            double share = 0.5 * edge_measure;
            node_area[n0] += share;
            node_area[n1] += share;

            // Acumular por elemento (SUMA, no max)
            int elem_id = faceList[i].elem_id;
            m_elem_area[elem_id] += edge_measure;
        }
    }
  }  
  else {
  ////////////////////////////////////
  //////// CALCULATE AREA (FOR CONTACT)
  // Allocate and initialize nodal area array

  // Compute nodal areas from external triangular faces
  for (int i = 0; i < m_faceCount; i++) {
      if (faceList[i].count == 1) { // External face
          int n0 = faceList[i].nodes[0];
          int n1 = faceList[i].nodes[1];
          int n2 = faceList[i].nodes[2];

          // Get coordinates of the nodes
          double3 p0 = getPosVec3(n0);
          double3 p1 = getPosVec3(n1);
          double3 p2 = getPosVec3(n2);

          // Compute face area via cross product
          double3 v1 = p1 - p0;
          double3 v2 = p2 - p0;
          double3 cross_ = cross(v1,v2);
          double area = 0.5 * norm(cross_);

          // Distribute area equally to the three nodes
          double area_share = area / 3.0;
          node_area[n0] += area_share;
          node_area[n1] += area_share;
          node_area[n2] += area_share;

          // if (abs(p0.z - 0.03)<5.e-4 || abs(p1.z - 0.03)<5.e-4 || abs(p2.z - 0.03)<5.e-4){
            // top_faces++;          
            // printf("FACE NORMAL %.4e %.4e %.4e\n",cross_.x,cross_.y,cross_.z);
          // }
          int elem_id = faceList[i].elem_id;
          if (!elem_flags[elem_id]){
            m_elem_area[elem_id] = area;
            elem_flags[elem_id] = true;
            //if (abs(p0.z - 0.03)<5.e-4 || abs(p1.z - 0.03)<5.e-4 || abs(p2.z - 0.03)<5.e-4)
            //top_area +=area;

          } 
          else {
            if (area > m_elem_area[elem_id]) m_elem_area[elem_id] = area;
          }
      }//==1 
    }//Face Count
    //printf ("Top Faces: %d\n",top_faces);
    //printf("--------------------------TOP AREA %.4e\n",top_area);
  }
  #ifdef CUDA_BUILD
    free_t(elem_flags);
  #endif
  
}

/////NEW
/////REFACTORING 
//~ void Domain_d::SearchExtNodes() {
    //~ // Initialize
    //~ Face *faceList = nullptr;
    //~ malloc_t(faceList, Face, m_elem_count*ELFAC);
    //~ int faceCount = 0;
    
    //~ // Process elements
    //~ for (int i = 0; i < m_elem_count; i++) {
        //~ int elements[ELNOD];
        //~ for(int ne=0; ne<m_nodxelem; ne++) {
            //~ elements[ne] = m_elnod[m_nodxelem*i+ne];
        //~ }
        //~ addTriangleFaces(faceList, faceCount, elements);
    //~ }

    //~ // Find external nodes
    //~ for(int i = 0; i < faceCount; i++) {
        //~ if(faceList[i].count == 1) {
            //~ for(int j = 0; j < FACENOD; j++) {
                //~ ext_nodes[faceList[i].nodes[j]] = true;
            //~ }
        //~ }
    //~ }

    //~ // Calculate areas with validation
    //~ for(int i = 0; i < faceCount; i++) {
        //~ if(faceList[i].count == 1) {
            //~ int n[3] = {faceList[i].nodes[0], faceList[i].nodes[1], faceList[i].nodes[2]};
            //~ double3 p[3] = {getPosVec3(n[0]), getPosVec3(n[1]), getPosVec3(n[2])};
            
            //~ double3 edge1 = p[1] - p[0];
            //~ double3 edge2 = p[2] - p[0];
            //~ double3 normal = cross(edge1, edge2);
            //~ double area = 0.5 * norm(normal);
            
            //~ if(area > 1.0) {
                //~ printf("Large area face %d: nodes (%d,%d,%d) area %e\n",
                      //~ i, n[0], n[1], n[2], area);
                //~ printf("Coordinates:\n");
                //~ for(int k=0; k<3; k++) 
                    //~ printf("  Node %d: (%f, %f, %f)\n", n[k], p[k].x, p[k].y, p[k].z);
            //~ }
            
            //~ double third_area = area / 3.0;
            //~ for(int k=0; k<3; k++)
                //~ node_area[n[k]] += third_area;
        //~ }
    //~ }
    
    //~ free_t(faceList);
//~ }

struct AssignValue {
    __device__ void operator()(double& x) const {
        x = 300.0e6;
    }
};

///// FUNCTOR SHOULD BE DEVICE
struct AssignVal {
    __device__ void operator()(double& x) const {
        x = 300.0e6;
        printf("VAL %.6e\n",x);
    }
};

dev_t void Domain_d::InitElemValues(double *arr, double val){
  
  par_loop(i, m_elem_count)
    arr[i] = val;
  
}

//THIS ASSUMES 
void Domain_d::InitValues(){


  
  #ifdef CUDA_BUILD
    //AssignValueFunctor assignValueFunctor(1.0e10);
    //parallel::for_each(sigma_y, sigma_y + m_elem_count, assignValueFunctor);


    //InitElemValuesKernel<<<1,1>>>(this,sigma_y, 300.0e6); //Another approach
    parallel::for_each(this->u, this->u + m_dim*m_node_count, AssignValueFunctor(0.0));
    parallel::for_each(this->v, this->v + m_dim*m_node_count, AssignValueFunctor(0.0));
    
    //InitElemValuesKernel<<<1,1>>>(this,this->v, 0.0);
    
    InitStressesFromMatKernel<<<1,1>>>(this);
  #else
    initElemArrayCPU(this,pl_strain,1,0.0)
    //initElemArrayCPU (this,sigma_y,1,1.0e10)  
    for (int e=0;e<m_elem_count;e++){
      sigma_y[e] = mat[e]->sy0;      
    }
  #endif
}

//Serial in order to be called as <<<1,1>>>
dev_t void Domain_d::InitStressesFromMat(){
    //par_loop(e,m_elem_count){
  //NOT PARALLEL
  for (int e=0;e<m_elem_count;e++)    
      sigma_y[e] = mat[e]->sy0;      
      
  //}
  
}


///// THIS SHOULD BE FOLLOWED BY ALLOCATEBCS()
///// TODO: ASSIGN ZONE FIRST //////
int Domain_d::AddBCVelZone(const vector_t &start, const vector_t &end, const vector_t &vel){
  int partcount = 0;
  for (size_t a=0; a<m_node_count; a++){
    bool included=true;
      if (x[m_dim*a] < start.x || x[m_dim*a] > end.x)
        included = false;
      if (x[m_dim*a+1] < start.y || x[m_dim*a+1] > end.y)
        included = false;
      if (m_dim > 2)
        if (x[m_dim*a+2] < start.z || x[m_dim*a+2] > end.z)
          included = false;
     
    if (included){     
      AddBCVelNode(a,0,vel.x);      AddBCVelNode(a,1,vel.y);      
      if (m_dim>2) AddBCVelNode(a,2,vel.z);
      partcount++;
    }
  }

  return partcount;
}


///// IMPORTANT: m_gp_count 

void Domain_d::SetDimension(const int &node_count, const int &elem_count){
  
  m_node_count = node_count;
  m_elem_count = elem_count;
  
  // NODAL VARIABLES
  cout << "Setting explicit dimensions for "<< node_count<<"nodes and "<< m_dim<<"dim"<<endl;
  malloc_t (x,      double,node_count*m_dim);
  malloc_t (v,      double,node_count*m_dim);
  malloc_t (a,      double,node_count*m_dim);
  malloc_t (u,      double,node_count*m_dim);
  malloc_t (u_dt,   double,node_count*m_dim);
  
  malloc_t (prev_a, double,node_count*m_dim);  
	//cudaMalloc((void **)&m_f, node_count * sizeof (double) * 3);
  malloc_t (m_fi,double,node_count*m_dim); //Internal forces
  malloc_t (m_fe,double,node_count*m_dim);

  //cout << "mdiag "<<endl;  
  malloc_t (m_mdiag, double,node_count);
  //malloc_t (m_mglob, double,node_count*node_count); //TODO: MAKE SPARSE. DEALLOCATED AFER DIAG CALCULATION

  //////if thermal

  malloc_t (T,                double, node_count);
  malloc_t(m_dTedt,           double, m_elem_count * m_nodxelem);
  malloc_t(ps_energy,         double, m_elem_count);
  malloc_t (q_cont_conv,      double, node_count);
  malloc_t (node_area,        double, node_count); /////USED FOIR THERMAL CONTACT
  //TODO: MAYBE THIS COULD BE DISCONTINUED 
  malloc_t (p_node,        double, node_count); /////USED FOIR THERMAL CONTACT

  
  //Used for contact ()
  malloc_t (m_elem_length,        	double, m_elem_count); /////USED FOIR THERMAL CONTACT
  malloc_t (m_mesh_in_contact,      int, 	m_node_count); /////USED FOIR CONTACT

  malloc_t (m_elem_min_angle,        	double, m_elem_count); /////USED FOIR THERMAL CONTACT
  malloc_t (m_elem_max_angle,        	double, m_elem_count); /////USED FOIR THERMAL CONTACT
  
  /// MATRICES ///
  /// dHxy_detJ: DIM X NODX ELEM
  /// TO AVOID EXCESSIVE OFFSET, SPLIT DIMENSIONS
  //cudaMalloc((void **)&m_dH_detJ_dx, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));
  //cudaMalloc((void **)&m_dH_detJ_dy, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));  
  //cudaMalloc((void **)&m_dH_detJ_dz, m_nodxelem * m_elem_count * m_gp_count * sizeof (double));  
  malloc_t (m_H,          double, m_dim * m_nodxelem * m_elem_count * m_gp_count);
  
  cout <<"setting deriv allocation size of "<<m_dim * m_nodxelem * m_elem_count * m_gp_count<<endl;
  cout <<"mdim"<<m_dim << " "<<m_nodxelem<<" "<< m_elem_count <<" "<<m_gp_count<<endl;
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
  malloc_t(m_eps,       double, 6 * m_elem_count * m_gp_count ); 

  malloc_t(m_str_rate_prev,  double, 6 * m_elem_count * m_gp_count );   
  
  //FOR PLASTIC STRAIN HEAT GENERATION
  malloc_t(m_strain_pl_incr,  double, 6 * m_elem_count * m_gp_count );     
  malloc_t(m_q_plheat,        double,     m_elem_count * m_gp_count );       
  
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

  malloc_t(ut_prev, double,node_count*m_dim);  
  
  
  m_contsurf_count = 0;

  //plastic
  malloc_t (pl_strain, double, m_elem_count * m_gp_count);
  malloc_t (sigma_y, double, m_elem_count * m_gp_count);  
  
  malloc_t (m_elem_area, double, m_elem_count);
  //MODIFY THIS!!! IS A LOT OF SPACE
  malloc_t(faceList, Face, m_elem_count*ELFAC);
  
  
  m_pl_energy = 0.0;
  
  malloc_t(m_elem_neigh,        int, m_elem_count*m_nodxelem);
  malloc_t(m_elem_neigh_count,  int, m_elem_count);
  
  
    /// IMPLICIT THINGS
  //malloc_t(m_Kmat, Matrix, m_elem_count );  Written for asinglepointer
  //~ m_Kmat = new Matrix*[m_elem_count];  // or use malloc_t macro if it's defined
  //~ for (int e=0;e<m_elem_count;e++)
    //~ m_Kmat[e] = new Matrix(m_nodxelem* m_dim,m_nodxelem* m_dim);


  //REMESHING
  malloc_t(m_sigma_prev,     double, 6 * m_elem_count * m_gp_count );   
  malloc_t(m_tau_prev,     double, 6 * m_elem_count * m_gp_count );   
  malloc_t(pl_strain_prev,   double, 1 * m_elem_count * m_gp_count );     
  malloc_t (m_vprev,          double,node_count*m_dim);

  
  if (m_dim == 2) local_faces = quad_edges; /////TODO: CHANGE TO DYNAMIC 
  
  if (m_dim==2){
  cout << " ALLOCATING HOURGLASS "<< m_elem_count*m_nodxelem*m_dim<<endl;
  malloc_t(m_hg_q, double, m_elem_count*m_nodxelem*m_dim); //m_hg_q[ m_elem_count ][ jmax ][ m_dim ]
  
  for (int i=0;i<m_elem_count*m_nodxelem*m_dim;i++)
    m_hg_q[i] = 0.0;
  }
  
  
}

void Domain_d::SetDimensionImplicit(const int &node_count, const int &elem_count){
  m_dim = 3;
  m_nodxelem = 4;
  m_gp_count = 1;
  cout << "Setting Dimension for implicit domain"<<endl;
  m_node_count = node_count;
  m_elem_count = elem_count;
  cout << "m_dim"<<endl;
  // NODAL VARIABLES
  
  malloc_t (x,      double,node_count*m_dim);
  malloc_t (v,      double,node_count*m_dim);
  malloc_t (a,      double,node_count*m_dim);
  malloc_t (u,      double,node_count*m_dim);
  malloc_t (u_dt,   double,node_count*m_dim);
  
  malloc_t (x_old,      double,node_count*m_dim);
  
  malloc_t (prev_a, double,node_count*m_dim);  
  malloc_t (prev_v, double,node_count*m_dim);  ///////ONLY USED ON IMPLICIT
  
	//cudaMalloc((void **)&m_f, node_count * sizeof (double) * 3);
  malloc_t (m_fi,double,node_count*m_dim); //Internal forces
  malloc_t (m_fe,double,node_count*m_dim);
  
  malloc_t (m_mdiag, double,node_count);
  //malloc_t (m_mglob, double,node_count*node_count); //TODO: MAKE SPARSE. DEALLOCATED AFER DIAG CALCULATION

  //////if thermal
  
  malloc_t (T,      double,node_count);
  malloc_t(m_dTedt,    double, m_elem_count * m_nodxelem);

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
  malloc_t(m_eps,       double, 6 * m_elem_count * m_gp_count );   
  //FOR PLASTIC STRAIN HEAT GENERATION
  malloc_t(m_strain_pl_incr,  double, 6 * m_elem_count * m_gp_count );     
  malloc_t(m_q_plheat,        double,     m_elem_count * m_gp_count );       
  
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
  malloc_t (node_area,        double, node_count); /////USED FOIR THERMAL CONTACT
  malloc_t(ut_prev, double,node_count*m_dim);  
  
  //TODO: MAYBE THIS COULD BE DISCONTINUED 
  malloc_t (p_node,        double, node_count); /////USED FOIR THERMAL CONTACT
  
  m_contsurf_count = 0;

  //plastic
  malloc_t (pl_strain, double, m_elem_count * m_gp_count);
  malloc_t (sigma_y, double, m_elem_count * m_gp_count);  
  
  /// IMPLICIT THINGS
  //malloc_t(m_Kmat, Matrix, m_elem_count );  Written for asinglepointer
  #ifndef CUDA_BUILD
  m_Kmat = new Matrix*[m_elem_count];  // or use malloc_t macro if it's defined
  for (int e=0;e<m_elem_count;e++)
    m_Kmat[e] = new Matrix(m_nodxelem* m_dim,m_nodxelem* m_dim);
  #else
  
  #endif
  
  malloc_t (m_elem_area, double, m_elem_count);
  //MODIFY THIS!!! IS A LOT OF SPACE
  malloc_t(faceList, Face, m_elem_count*ELFAC);
  
  malloc_t(m_str_rate_prev,  double, 6 * m_elem_count * m_gp_count );   
   
  
  m_pl_energy = 0.0;
  
  malloc_t(m_elem_neigh,        int, m_elem_count*4);
  malloc_t(m_elem_neigh_count,  int, m_elem_count);
  
  ////// CONTACT
  malloc_t (q_cont_conv,      double, node_count);
  malloc_t (node_area,        double, node_count); /////USED FOIR THERMAL CONTACT

  //Used for contact ()
  malloc_t (m_elem_length,        	double, m_elem_count); /////USED FOIR THERMAL CONTACT
  malloc_t (m_mesh_in_contact,      int, 	m_node_count); /////USED FOIR CONTACT


  //REMESHING
  malloc_t(m_sigma_prev,     double, 6 * m_elem_count * m_gp_count );   
  malloc_t(m_tau_prev,     double, 6 * m_elem_count * m_gp_count );   
  malloc_t(pl_strain_prev,   double, 1 * m_elem_count * m_gp_count );     
  malloc_t (m_vprev,      double,node_count*m_dim);
}

void Domain_d::SetDimensionExplicit(const int &node_count, const int &elem_count){
  
  
}

void Domain_d::Free(){
  free_t (x);
  free_t (v);
  free_t (a);
  free_t (u);
  free_t (u_dt);
  
  free_t (prev_a);  
	//cudaMalloc((void **)&m_f, node_count * sizeof (double) * 3);
  
  free_t (m_fi); //Internal forces
  free_t (m_fe);
  
  free_t (m_mdiag);
  //free_t (m_mglob); //TODO: MAKE SPARSE. DEALLOCATED AFER DIAG CALCULATION

  free_t (T);
  free_t(m_dTedt);
  free_t(m_q_plheat);

 
  free_t (m_H);

  free_t (m_dH_detJ_dx);
  free_t (m_dH_detJ_dy);
  free_t (m_dH_detJ_dz);


  //DOUBLE POINTERS ATTENTION
  // free_t (m_dHrs);     //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  // free_t (x2);         //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  // free_t (dH_dxyz); 
  
  free_t(elnod_h); ////// USED TO COMPUTE GLOBAL M MATRIX WHICH IS COMPUTED IN CPU (TODO: CHANGE)       
  free_t(dHxy_detJ ); ////NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX

  

  free_t(m_detJ );    
  
  free_t(m_nodxelem_e);
  
  free_t(m_ematm); //Elemental mas matrices

  free_t(m_str_rate);   
  free_t(m_rot_rate );     
  free_t(m_sigma);   
  free_t(m_tau); 
  free_t(m_eps);    
 

  free_t(p); 
  free_t(rho);   
  free_t(rho_0); 
  
  free_t(vol); 
  free_t(vol_0); 

  free_t(m_voln); 
  free_t(m_vol_0n); 

  free_t(m_f_elem);   
  free_t(m_f_elem_hg);   

  free_t(mat);   

  // AXISYMM
  free_t (m_radius);

  free_t(m_jacob);
    

  free_t(ext_nodes);
  free_t(contforce);

  free_t (pl_strain);
  free_t (sigma_y);  

  free_t(m_nodel);           // NODE ELEMENTS: ELEMENTS SHARED BY EACH NODE [nodxelem* node_count] call:  m_nodel[nod,elcount] =EL i,e m_nodel[nod_offset+elcount]
  free_t(m_nodel_loc);        //
  free_t(m_nodel_offset);    //OFFSET OF THE
  free_t(m_nodel_count);    
  free_t(m_elnod_count);   /// FOR CONTACT, TO REPLACE FOR m_node_count
  
  // unsigned int    *m_contsurf_count;
  // unsigned int    *m_contsurf_elemcount;   //FOR EACH OF THE ABOVE  
  // unsigned int    *m_contsurf_elem;        //ELEMENT POS OF THE CONTACT ELEMENT 

  // ////////////////////// CONTACT 
	// // TODO, EACH RIGID PARTICLE SHOULD 
  // int   *contelem; //ELEMENT OF TRIMESH FROM "RIGID" PARTICLE, ALL FIRST PARTICLES ARE ZERO
  // TriMesh_d *trimesh;
  // int trimesh_count;
  // int *mesh_id; //particle mesh ID	
	// bool *ext_nodes;
  // int ext_nodes_count;
  // double *contforce; 
  // bool contact;
  
  //free_t(bcx_nod);  free_t(bcy_nod);  free_t(bcz_nod);
  //free_t(bcx_val);  free_t(bcy_val);  free_t(bcz_val);
  //bc_count[0] = bc_count[1] = bc_count[2] = 0;

  free_t(m_elem_area);
  free_t(faceList);    
  free_t(node_area);
  
  
  free_t (m_elem_neigh);
  free_t (m_elem_neigh_count);  
  
  free_t(m_elem_min_angle);
  free_t(m_elem_max_angle);
  free_t (m_mesh_in_contact);

}


void Domain_d::AssignMaterial (Material_ *material_h) {
//    cudaMalloc((void**)&materials, 1 * sizeof(Material_ )); //
    malloc_t(materials, Material_,1);
    memcpy_t(materials, material_h, 1 * sizeof(Material_));	
}
#ifdef CUDA_BUILD

struct Increment {
    __device__ void operator()(double& x) {
        x += 1.0;
    }
};

#else
struct Increment {
    void operator()(double& x) {
        x += 1.0;
    }
};

#endif

//// NEW PARALLEL INCREMENT
void Domain_d::AssignMatAddress_(){
    int N = 100;
    double* d_array;
/*
    // Allocate device memory
    malloc_t(&d_array, N * sizeof(double));
    
    // Initialize with zeros
    cudaMemset(d_array, 0, N * sizeof(double));
  parallel::for_each(d_array, d_array + N, Increment());

    // Free memory
    cudaFree(d_array);

    std::cout << "CUDA for_each completed!" << std::endl;
    //return 0;
    */
}

dev_t void Domain_d::AssignMatAddress(){
  par_loop(i, m_elem_count)
    mat[i] = &materials[0];
  
}

///// MAKE TEMPLATE SET VECTOR
void Domain_d::setDensity(const double &r){
  double *rho_h = new double [m_elem_count];
  for (int n=0;n<m_elem_count;n++) rho_h[n] = r;
  memcpy_t(this->rho_0, rho_h, sizeof(double) * m_elem_count);    
  
  delete rho_h;
  
}

dev_t void Domain_d::UpdatePrediction(){
  par_loop (i,m_node_count){
      for (int j = 0; j < m_dim; j++) {
          //u_[i][j] = dt * (v_[i][j] + (0.5 - m_beta) * dt * prev_a_[i][j]);
          //NEW; global
          int ig = i*m_dim + j; //BY NOW is 2D
          
          u_dt[ig] = dt * (v[ig] + (0.5 - m_beta) * dt * prev_a[ig]);
          v   [ig] += (1.0 - m_gamma) * dt * prev_a[ig];    

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
  
    //for (int i = 0; i < m_node_count; i++) {
    par_loop (i,m_node_count) {
        for (int j = 0; j < m_dim; j++) {
            int ig = i*m_dim + j;

            //printf ("a ig %f, node %d, dim %d \n",  a[ig]);
            a[ig] = f * (a[ig]  -m_alpha * prev_a[ig]); //GLOBAL
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

//      for (int i = 0; i < m_node_count; i++) {
  par_loop (i,m_node_count) {
    for (int j = 0; j < m_dim; j++) {
        // u_[i][j] += m_beta * dt * dt * a_[i][j];
        // x_[i][j] += u_[i][j];

        int ig = i*m_dim + j;

        //printf("GLOBAL IND %d\n", ig);
        u_dt[ig] += m_beta * dt * dt * a[ig];
        x[ig] += u_dt[ig];
        prev_a[m_dim*i+j] = a[m_dim*i+j];
        u[m_dim*i+j] += u_dt[m_dim*i+j];         
    }
    //printf ("U Node %d %.6e %.6e %.6e\n", i, u[m_dim*i],u[m_dim*i+1],u[m_dim*i+2] );     
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
  //for (int d=0;d<m_dim;d++){cout << "Allocated "<<bc_count[d]<< " Velocity BCs for dof "<<d<<endl; }
  cout <<"Allocated BCS: "<< endl;
  cout << "X "<<bcx_nod_h.size()<<endl;  cout << "Y "<<bcy_nod_h.size()<<endl;  cout << "Z "<<bcz_nod_h.size()<<endl;
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
  //printf("Dim %d, BC Count %d\n", dim, bc_count[dim]);
  
  par_loop (n,bc_count[dim]){
    double val;
    //printf("thread %d, Imposing Vel in dim %d, %d Conditions, val %f\n", n, dim, bc_count[dim], bcx_val[n]);
    //printf("BCV dim %d\n", dim);
    // printf("VEL BC \n");
    if (dim == 0)       {/*printf ("dim %d node %d, val %f\n",dim,bcx_nod[n],bcx_val[n]);*/  v[m_dim*bcx_nod[n]+dim] = bcx_val[n]; }
    else if (dim == 1)  {/*printf ("dim %d node %d val %f \n",dim,bcy_nod[n],bcy_val[n]);*/  v[m_dim*bcy_nod[n]+dim] = bcy_val[n];}
    else if (dim == 2)  {/*printf ("dim %d node %d, val %f\n",dim,bcz_nod[n],bcz_val[n]);*/  v[m_dim*bcz_nod[n]+dim] = bcz_val[n]; }
  }
  
}

dev_t void Domain_d::ImposeBCA(const int dim){
  par_loop (n,bc_count[dim]){
    double val;
    //printf("thread %d, Imposing Vel in dim %d, %d Conditions\n", n, dim, bc_count[dim]);
    //printf("ACEL BC \n");
    if (dim == 0)       {/*printf ("dim %d val %f, Nod %d\n",dim, bcx_val[n],bcx_nod[n]);*/ a[m_dim*bcx_nod[n]+dim] = 0.0; }
    else if (dim == 1)  {/*printf ("dim %d val %f, Nod %d\n",dim, bcy_val[n],bcy_nod[n]); */a[m_dim*bcy_nod[n]+dim] = 0.0;}
    else if (dim == 2)  {/*printf ("dim %d val %f, Nod %d\n",dim, bcz_val[n],bcz_nod[n]); */a[m_dim*bcz_nod[n]+dim] = 0.0; }
  }
  
}

void Domain_d::AddBoxLength(vector_t const & V, vector_t const & L, const double &r,const bool &red_int, const bool &tritetra){
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
    if (!tritetra)
      m_nodxelem = 4;
    else 
      m_nodxelem = 3;
    cout << "ne "<<m_nodxelem<<endl;
    if (!red_int) m_gp_count = 4;
  } else {
    nel[2] = (int)(L.z/(2.0*r));
    nel[2] = (int)(L.z/(2.0*r));
    
    if (!tritetra)
      m_nodxelem = 8;
    else
      m_nodxelem = 4;
      
    if (!red_int) m_gp_count = 8; 
  }
  

  Xp.z = V.z ;
    

    // write (*,*) "Creating Mesh ...", "Elements ", neL.y, ", ",neL.z
  int nc;
  if (m_dim == 2) nc = (nel[0] +1) * (nel[1]+1);
  else            nc = (nel[0] +1) * (nel[1]+1) * (nel[2]+1);
  int ne;
  
  if (!tritetra)
    ne= nel[0]*nel[1]*nel[2];
  else {
    if (m_dim == 2)
      ne= 2 * nel[0]*nel[1]*nel[2];
    else
      ne= 6 * nel[0]*nel[1]*nel[2];
  }
  
  cout << "Mesh created. Element count: "<< nel[0]<<", "<<nel[1]<<", "<<nel[2]<<endl;
  cout << "Element nodes size "<<ne<<endl;
  //thisAllocateNodes((nel[0] +1) * (nel[1]+1) * (nel[2]+1));
    // print *, "Element count in XYZ: ", nel(:)
    // write (*,*) "Box Node count ", node_count

	this->SetDimension(nc,ne);	 //AFTER CREATING DOMAIN
  cout << "Mesh generated. Node count: " << nc<<". Element count: "<<ne<<endl;
  cout << "Dimension is: "<<m_dim<<endl;
  //SPH::Domain	dom;
	//vector_t *x =  (vector_t *)malloc(dom.Particles.size());

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
        x_h[m_dim*p  ] = Xp.x;  
        x_h[m_dim*p+1] = Xp.y;
        if (m_dim == 3) x_h[m_dim*p+2] = Xp.z;
        //nod%x(p,:) = Xp(:);
        //cout << "node " << p <<"X: "<<Xp.x<<"Y: "<<Xp.y<<"Z: "<<Xp.z<<endl;
        p++;
        Xp.x = Xp.x + 2.0 * r;
      }
      Xp.y = Xp.y + 2.0 * r;
    }// 
    Xp.z = Xp.z + 2 * r;

  //cout <<"m_node size"<<m_node.size()<<endl;
  } 
		memcpy_t(this->x,   x_h, m_dim*sizeof(double) * m_node_count);    
    printf("X\n");
    //printVec(this->x);
    //printf("x_h\n");
    //printVec(x_h);
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
    
    if (m_dim == 2) {

      if (!tritetra){
      int ei = 0;
      for (int ey = 0; ey < nel[1];ey++){
        for (int ex = 0; ex < nel[0];ex++){

        int nb1 = (nel[0]+1)* ey    + ex;    
        int nb2 = (nel[0]+1)*(ey+1) + ex;
        elnod_h[ei  ] = nb1;                        nodel_count_h[nb1  ] ++;             
        elnod_h[ei+1] = nb1 + 1;                    nodel_count_h[nb1+1] ++;             
        elnod_h[ei+2] = nb2 + 1;                    nodel_count_h[nb2+1] ++;      
        elnod_h[ei+3] = (nel[0]+1)*(ey+1) + ex;     nodel_count_h[nb2 ] ++;      
			
				//for (int i=0;i<m_nodxelem;i++)cout << elnod_h[ei+i]<<", ";
					//cout << "Nel x : "<<nel[0]<<endl;
					//cout << "nodes "<<endl;
					ei += m_nodxelem;
					 }
      }
      }else {//tritetra
        int ei = 0;
        for (int ey = 0; ey < nel[1];ey++){
          for (int ex = 0; ex < nel[0];ex++){

          int nb1 = (nel[0]+1)* ey    + ex;    
          int nb2 = (nel[0]+1)*(ey+1) + ex;
          // |\
          // | \
          // ----
          elnod_h[ei  ] = nb1;                        nodel_count_h[nb1  ] ++;             
          elnod_h[ei+1] = nb1 + 1;                    nodel_count_h[nb1+1] ++;             
          elnod_h[ei+2] = nb2;                        nodel_count_h[nb2]   ++;  
          ei += 3;
            
          elnod_h[ei  ] = nb1 + 1;                    nodel_count_h[nb1+1] ++;             
          elnod_h[ei+1] = nb2 + 1;                    nodel_count_h[nb2+1] ++;            
          elnod_h[ei+2] = nb2;                        nodel_count_h[nb2  ] ++;      
        
          //for (int i=0;i<m_nodxelem;i++)cout << elnod_h[ei+i]<<", ";
            //cout << "Nel x : "<<nel[0]<<endl;
            //cout << "nodes "<<endl;
          ei += 3; 
             }
        }      
      }//tritetra
    } else { //dim: 3
      if (!tritetra) {
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

             //cout << "Nel x : "<<nel[0]<<", Element: " << ei/m_nodxelem<<endl;
             //cout << "nodes "<<endl;
             //cout << "nodxelem " <<m_nodxelem<<endl;
             //cout << "nb1, nb2 "<< nb1 << ", "<<nb2<<" nnodz"<< nnodz*(ez+1)<<"ez "<<ez<<endl;
             for (int i=0;i<m_nodxelem;i++)cout << elnod_h[ei+i]<<", ";
             ei += m_nodxelem;

             }
        } 
      } else {
        int ei = 0; //ELEMENT INTERNAL NODE (GLOBAL INDEX)
        int nnodz = (nel[0]+1)*(nel[1]+1);
        
        //////https://www.researchgate.net/publication/329270810_There_are_174_subdivisions_of_the_hexahedron_into_tetrahedra
        //////FIGURE 9
        
        for (int ez = 0; ez < nel[2];ez++){
          for (int ey = 0; ey < nel[1];ey++){
            for (int ex = 0; ex < nel[0];ex++){
            
              int nb1 = nnodz*ez + (nel[0]+1) *  ey    + ex;
              int nb2 = nnodz*ez + (nel[0]+1) * (ey+1) + ex;   
              int nhex[] = {nb1,         nb1+1,        nb2+1, nb2, 
                                          nb1 + nnodz, nb1 + nnodz+1,nb2 + nnodz + 1, nb2 + nnodz};
              //cout << "HEXA NODES "<<endl;
              //for (int i=0;i<8;i++)
              //  cout << nhex[i]<<" ";
              //cout <<endl;
              
              //1st valid decomp
              //0,1,3,4
              //1,3,4,5
              //1,2,3,5,
              //3,4,5,7
              //2,3,5,6
              
              //Another valid decomp
              //0,1,3,7
             // 0,1,7,5
              //0,3,7,4
              //1,7,5,6
              //1,3,7,6
              //3,7,4,6
              

              elnod_h[ei] = nhex[0];elnod_h[ei+1] = nhex[1];elnod_h[ei+2] = nhex[3];elnod_h[ei+3] = nhex[5]; ei += m_nodxelem;
              elnod_h[ei] = nhex[1];elnod_h[ei+1] = nhex[2];elnod_h[ei+2] = nhex[3];elnod_h[ei+3] = nhex[5]; ei += m_nodxelem;
              elnod_h[ei] = nhex[0];elnod_h[ei+1] = nhex[5];elnod_h[ei+2] = nhex[3];elnod_h[ei+3] = nhex[4]; ei += m_nodxelem;
              elnod_h[ei] = nhex[4];elnod_h[ei+1] = nhex[5];elnod_h[ei+2] = nhex[3];elnod_h[ei+3] = nhex[7]; ei += m_nodxelem;
              elnod_h[ei] = nhex[5];elnod_h[ei+1] = nhex[6];elnod_h[ei+2] = nhex[3];elnod_h[ei+3] = nhex[7]; ei += m_nodxelem;
              elnod_h[ei] = nhex[5];elnod_h[ei+1] = nhex[2];elnod_h[ei+2] = nhex[3];elnod_h[ei+3] = nhex[6]; ei += m_nodxelem;
              
              nodel_count_h[nhex[0]] +=2;nodel_count_h[nhex[1]] +=2;   nodel_count_h[nhex[2]] +=2; nodel_count_h[nhex[3]] +=6;
              nodel_count_h[nhex[4]] +=2;nodel_count_h[nhex[5]] +=6;   nodel_count_h[nhex[6]] +=2; nodel_count_h[nhex[7]] +=2;                
              
              
              
            }
          }
        }
        
      }/////TRITETRA
		}//if dim 
    
    cout << "Done." <<endl;
    
    cout << "Allocating "<< m_elem_count<< " and "<<m_nodxelem<< "nodes x elem" <<endl;
    ///// FOR DIFFERENT ELMENT NODE COUNT
    int *m_nodxelem_eh  = new int [m_elem_count];
    int *m_elnodoffset_h = new int [m_elem_count];
    
    /*
    m_elnodoffset[0]=0;
    for (int e = 0;e<m_elem_count;e++){
      m_nodxelem_eh[e] = 4;
      m_elnodoffset_h[e] += m_nodxelem_eh[e];
    }
		memcpy_t(this->m_nodxelem_e, m_nodxelem_eh, sizeof(int) * m_elem_count); 
		memcpy_t(this->m_elnodoffset, m_elnodoffset, sizeof(int) * m_elem_count); 
    */
    
        
    
    //cudaMalloc((void **)&m_elnod, m_elem_count * m_nodxelem * sizeof (int));	
    malloc_t(m_elnod, unsigned int, m_elem_count * m_nodxelem);
    cout << "COPYING "<<m_elem_count * m_nodxelem<< " element nodes "<<endl;
		memcpy_t(this->m_elnod, elnod_h, sizeof(int) * m_elem_count * m_nodxelem); 
    
    // printf("ELNOD \n");
    // for (int e=0;e<m_elem_count;e++){
      // for (int n=0;n<m_nodxelem;n++)
        // printf("%d ",m_elnod[e*m_nodxelem+n]);
      // printf("\n");
    // }
    //cudaMalloc(&m_jacob,m_elem_count * sizeof(Matrix ));
    cout << "Done"<<endl;
    
    //////////////////// ELEMENT SHARED BY NODES (FOR PARALLEL NODAL MODE ASSEMBLY) ///////////////////////////////
    int nodel_tot = 0;
    for (int n=0;n<m_node_count;n++){
      nodel_offset_h[n] = nodel_tot;
      nodel_tot        += nodel_count_h[n];
      //cout << "NodEL tot " << nodel_tot<<endl;
      //cout << "Node "<< n << " Shared elements: "<<nodel_count_h[n]<<endl;

    }
    cout << "Size of Nodal shared Elements vector "<< nodel_tot<<endl;
		int *nodel_h       = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT
    int *nodel_loc_h   = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT    
    
    //Reset nodelcount (is incremented latere)
    for (int n=0;n<m_node_count;n++)  nodel_count_h[n] = 0;    
    
    //ALLOCATE NODEL, WHICH ARE ELEMEENTS SHARED BY A NODE
    //THIS IS FOR MORE EFFICIENT ELEMENTFORCES  ASSEMBLY
    //
    for (int e=0;e<m_elem_count;e++){
      int offset = m_nodxelem * e;
      for (int ne=0;ne<m_nodxelem;ne++){
        int n = elnod_h[offset+ne];
        if (nodel_offset_h[n] + nodel_count_h[n] > nodel_tot)
          cout << "ERRROR in node index,index "<<nodel_offset_h[n] + nodel_count_h[n]<<", node "<<n<<"element "<<e<<endl;
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
		  
    /*
    // ///// TESTING
    cout << endl<<endl;
    for (int n=0;n<m_node_count;n++){
      cout << "M node offset:"<<nodel_offset_h[n]<<endl;
      cout << "Node  "<< n << " Elements ("<< nodel_count_h[n]<<")"<<endl;
      for (int ne=0;ne<nodel_count_h[n];ne++) cout << nodel_h[nodel_offset_h[n]+ne]<<", ";
      cout << endl;
      cout << "Node  "<< n << " Elements Internal Node"<<endl;
      for (int ne=0;ne<nodel_count_h[n];ne++) cout << nodel_loc_h[nodel_offset_h[n]+ne]<<", ";
      cout << endl;
    }
    */
    
    cout << "Mesh generation done. "<<endl;
		m_tot_mass = L.x*L.y;
    if (m_dim == 3){
      m_tot_mass*=L.z;
    }

		delete [] /*elnod_h, */nodel_count_h, nodel_h, nodel_loc_h,nodel_offset_h;
}


//USED FOR PARALLEL ASSEMBLY AND AVERAGING
void Domain_d::setNodElem(int *elnod_h){


     int *nodel_count_h  = new int [m_node_count];
     int *nodel_offset_h = new int [m_node_count];
    for (int n=0;n<m_node_count;n++){
      nodel_count_h[n] = 0;
    }
    
    //Must initialize nodel_count_h
  int offset = 0;
  for (int e=0;e<m_elem_count;e++){
    //cout << "Element "<<e<<endl;
    for (int ne=0;ne<m_nodxelem;ne++){
      if (elnod_h[offset+ne]<m_node_count){
        //if (elnod_h[offset+ne]==3)
          //cout << "Node 3 shared element "<< e<<endl;
        nodel_count_h[elnod_h[offset+ne]]++;
      }else 
        cout << "ERROR setting node element, element "<<e <<", node "<<ne<<", global "<<elnod_h[offset+ne]<<"> Node Count "<<endl;
    }
    offset+=m_nodxelem;
  }
  cout << "Allocating "<< m_elem_count<< " and "<<m_nodxelem<< "nodes x elem" <<endl;
  ///// FOR DIFFERENT ELMENT NODE COUNT
  int *m_nodxelem_eh  = new int [m_elem_count];
  int *m_elnodoffset_h = new int [m_elem_count];
  

  
  //cudaMalloc((void **)&m_elnod, m_elem_count * m_nodxelem * sizeof (int));	
  malloc_t(m_elnod, unsigned int, m_elem_count * m_nodxelem);
  cout << "COPYING "<<m_elem_count * m_nodxelem<< " element nodes "<<endl;
  memcpy_t(this->m_elnod, elnod_h, sizeof(int) * m_elem_count * m_nodxelem); 
  cout << "Done"<<endl;
  
  //////////////////// ELEMENT SHARED BY NODES (FOR PARALLEL NODAL MODE ASSEMBLY) ///////////////////////////////
  int nodel_tot = 0;
  for (int n=0;n<m_node_count;n++){
    nodel_offset_h[n] = nodel_tot;
    nodel_tot        += nodel_count_h[n];
    //cout << "NodEL tot " << nodel_tot<<endl;
    //cout << "Node "<< n << " Shared elements: "<<nodel_count_h[n]<<endl;

  }
  cout << "Size of Nodal shared Elements vector "<< nodel_tot<<endl;
  int *nodel_h       = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT
  int *nodel_loc_h   = new int [nodel_tot];          //ASSUMED EACH NODE SHARES 8 ELEMENT    
  
  //Reset nodelcount (is incremented latere)
  for (int n=0;n<m_node_count;n++)  nodel_count_h[n] = 0;    
  
  //ALLOCATE NODEL, WHICH ARE ELEMEENTS SHARED BY A NODE
  //THIS IS FOR MORE EFFICIENT ELEMENTFORCES  ASSEMBLY
  //
  for (int e=0;e<m_elem_count;e++){
    int offset = m_nodxelem * e;
    for (int ne=0;ne<m_nodxelem;ne++){
      int n = elnod_h[offset+ne];
      if (nodel_offset_h[n] + nodel_count_h[n] > nodel_tot)
        cout << "ERRROR in node index,index "<<nodel_offset_h[n] + nodel_count_h[n]<<", node "<<n<<"element "<<e<<endl;
      nodel_h     [nodel_offset_h[n] + nodel_count_h[n]] = e;
      nodel_loc_h [nodel_offset_h[n] + nodel_count_h[n]] = ne;
      
      nodel_count_h[n]++;
    }//nod x elem 
  }
  
  malloc_t (m_nodel,        int,nodel_tot);
  malloc_t (m_nodel_loc,    int,nodel_tot);
  
  malloc_t (m_nodel_offset, int, m_node_count);
  malloc_t (m_nodel_count,  int, m_node_count);
  
  //THIS IS ONLY FOR COMPLETENESS IN CASE OF CPU, SHOULD BE BETTER TO WRITE ON FINALL ARRAY
  memcpy_t(this->m_nodel,         nodel_h,        sizeof(int) * nodel_tot); 
  memcpy_t(this->m_nodel_loc,     nodel_loc_h,    sizeof(int) * nodel_tot);
   
  memcpy_t(this->m_nodel_offset,  nodel_offset_h, sizeof(int) * m_node_count);  //OFFSET FOR PREVIOUS ARRAYS
  memcpy_t(this->m_nodel_count,    nodel_count_h, sizeof(int) * m_node_count);  //OFFSET FOR PREVIOUS ARRAYS
    
  /*
  // ///// TESTING
  cout << endl<<endl;
  for (int n=0;n<m_node_count;n++){
    cout << "M node offset:"<<nodel_offset_h[n]<<endl;
    cout << "Node  "<< n << " Elements ("<< nodel_count_h[n]<<")"<<endl;
    for (int ne=0;ne<nodel_count_h[n];ne++) cout << nodel_h[nodel_offset_h[n]+ne]<<", ";
    cout << endl;
    cout << "Node  "<< n << " Elements Internal Node"<<endl;
    for (int ne=0;ne<nodel_count_h[n];ne++) cout << nodel_loc_h[nodel_offset_h[n]+ne]<<", ";
    cout << endl;
  }
  */
  
  cout << "Node Elemets set. "<<endl;

    
  delete [] nodel_count_h;
  delete [] nodel_h;
  delete [] nodel_loc_h;
  delete [] nodel_offset_h;

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


void Domain_d::CreateFromLSDyna(lsdynaReader &reader){
 
  m_dim = 3;
  vector_t Xp;
  ///IMPORTANT BEFORE SETDIMENSION, m_gp_count and m_nodxelem must be set
  m_gp_count = 1;
  m_nodxelem = reader.m_elem[0].node.size(); //TO CHANGE (CONSTANT)  
  if (m_timeint_type == TimeInt::IMPLICIT)
    this->SetDimensionImplicit(reader.m_node.size(),reader.m_elem_count);	 //AFTER CREATING DOMAIN
  else 
    this->SetDimension(reader.m_node.size(),reader.m_elem_count);	 //AFTER CREATING DOMAIN
  
  cout << "Allocating "<< reader.m_node.size()<< " nodes "<<endl;
  double *x_H =  new double [m_dim*m_node_count];
  for (int n=0;n<m_node_count;n++){
    //cout << "Node "<<n<<endl;
    for (int d=0;d<3;d++){
      //cout <<reader.m_node[n].m_x[d]<< " ";
      x_H[3*n+d] = reader.m_node[n].m_x[d]; 
    }
    //cout <<endl;
  }
  memcpy_t(this->x,   x_H, m_dim*sizeof(double) * m_node_count);    
  delete[] x_H;
  

  int *elnod_h       = new int [m_elem_count * m_nodxelem]; //Flattened  
  cout << "Allocating "<< m_elem_count <<" elements, node x element: "<<m_nodxelem<<endl; //TO CHANGE (CONSTANT)
  
  int offs = 0;
  for (int e=0;e<reader.m_elem_count;e++){
    //cout << "Element node count "<<reader.m_elem[e].node.size()<<endl;
    for (int en=0;en<reader.m_elem[e].node.size();en++  ){
      elnod_h [offs+en] = reader.m_elem[e].node[en];
      //ls_node n = reader.getElemNode(e, en);
      //cout << reader.m_elem[e].node[en]<<" ";
      //cout << "Node id "<<n.m_id<<", xyz:"<<n.m_x[0]<<", "<<n.m_x[1]<<", "<<n.m_x[2]<<endl;
    }
    //cout << endl;
    offs += reader.m_elem[e].node.size();    
  }

  cout << "Element nodes "<<elnod_h [0]<<" "<<elnod_h [1] << " "<<elnod_h [2]<<endl;
    
  cout << "Setting elements per node "<<endl; 
  setNodElem(elnod_h);
  
  cout << "Node Size: "<<m_node_count<<endl;  
  cout << "Element Size: "<<m_elem_count<<endl;  
  
  delete [] elnod_h;
  
}

dev_t void Domain_d::calcElemJAndDerivatives () {
	////printf("calculating\n");
	////printf ("threadIdx.x %d, blockDim.x%d, blockIdx.x %d\n",threadIdx.x ,blockDim.x , blockIdx.x);
  //printf("allocating m_dim %d, n_nodxelem %d\n", m_dim, m_nodxelem);
  

   
  //printf("done\n");
   
	////printf ("e %d, elem_count %d\n",m_elem_count);
  par_loop (e, m_elem_count) {

  Matrix jacob(m_dim, m_dim);
  Matrix inv_j(m_dim, m_dim);
  Matrix x2(m_nodxelem, m_dim);  
  Matrix dHxy_detJ_loc(m_dim, m_nodxelem);
  
	//Matrix *jacob = new Matrix(m_dim, m_dim);    
	//Matrix *inv_j = new Matrix(m_dim, m_dim);    
  //Matrix *dHrs = new Matrix(m_dim, m_nodxelem);   /////////////////////////////// IF CREATION IS DYNAMIC ! (TEST IF )
  //Matrix *x2 = new Matrix(m_nodxelem, m_dim);  
  //Matrix *dHxy_detJ_loc = new Matrix(m_dim, m_nodxelem);

  int offset = m_gp_count * e;

  // integer :: e
  // ! !rg=gauss[ig]
  // ! !sg=gauss[jg]
  // real(fp_kind), dimension(dim,m_nodxelem) :: dHrs !!! USED ONLY FOR SEVERAL GAUSS POINTS
	//printf ("m_dim %d, nod x elem %d", m_dim, m_nodxelem);

  
  //cudaMalloc((void**)&dHrs_p, sizeof(Matrix));
	////printf("test %lf",dHrs.m_data[0]);
	//double dHrs_fl[m_dim* m_nodxelem];
	//dHrs.Print();

  //printf("x2 dimensions %d X %d\n", m_nodxelem, m_dim);
   
   // //printf("Jacob\n");jacob.Print();
   double gpc[8][3];

	//printf ("Matrices created\n");

      // do i=1,nodxelem
          // !print *, "elnod " , elem%elnod(e,i)
          // x2(i,:)=nod%x(elem%elnod(e,i),:)
      // end do
  int nind = e * m_nodxelem;
  for (int i=0;i<m_nodxelem;i++){
      ////TEMPLATIZE
      if (m_dim == 2){
        double2 x_ = Ptr_vector2(x,m_elnod[nind+i]);
        x2.Set(i,0,x_.x); x2.Set(i,1,x_.y); 
        //printf("x2\n");
        //x2.Print();
      } else {
        vector_t x_ = Ptr_vector_t(x,m_elnod[nind+i]); 
        x2.Set(i,0,x_.x); x2.Set(i,1,x_.y);        
        x2.Set(i,2,x_.z);
       
      }

        
      
      ////printf ("elnod %d, %lf %lf %lf \n",m_elnod[nind+i],x[m_elnod[nind+i]].x,x[m_elnod[nind+i]].y,x[m_elnod[nind+i]].z);
  } 
  //printf("x2\n");x2.Print();
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
        if (m_nodxelem == 4){
          for (int d=0;d<2;d++){
            // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)
            // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)
            // elem%jacob(e,gp,:,:) = 0.25*elem%jacob(e,gp,:,:)
            jacob.Set(0,d,0.25*(-x2.getVal(0,d)+x2.getVal(1,d)+x2.getVal(2,d)-x2.getVal(3,d))); 
            jacob.Set(1,d,0.25*(-x2.getVal(0,d)-x2.getVal(1,d)+x2.getVal(2,d)+x2.getVal(3,d)));
          }
        
          AdjMat(jacob, &inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
          //printf(" J ptr\n");
          //jacob.Print();
          //printf("ADJ J ptr\n");
          //inv_j.Print();          //printf("jacob\n");jacob.Print();
          //invj x dHdrs [-1,1,1,-1;  -1.-1,1,1] 
          for (int d=0;d<2;d++){        
            dHxy_detJ_loc.Set(d,0,0.25*(-inv_j.getVal(d,0)-inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,1,0.25*(inv_j.getVal(d,0)-inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,2,0.25*( inv_j.getVal(d,0)+inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,3,0.25*(-inv_j.getVal(d,0)+inv_j.getVal(d,1)));     
          }
        //dHxy_detJ_loc.Mul(0.25);
        } else if (m_nodxelem == 3){ //TRIANGLE CONSTANT ELEMENT
          //BENSON 2.4.5.2 N1 = r , N2 = s, n3 = 1 - r -s
           //dHdrs [1,0,-1;  -0,1,-1] x X2
          for (int d=0;d<2;d++){
            jacob.Set(0,d,(x2.getVal(0,d)-x2.getVal(2,d))); 
            jacob.Set(1,d,(x2.getVal(1,d)-x2.getVal(2,d)));          
          }
          AdjMat(jacob, &inv_j);
          //invj x dHdrs [-1, 1, 0,-1;  
          //              -1, 0, 1, 0]
          //                  
          for (int d=0;d<2;d++){    //col of dHdrs     
            dHxy_detJ_loc.Set(d,0,(inv_j.getVal(d,0)));   //row 1 of jacobian  
            dHxy_detJ_loc.Set(d,1,(inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,2,(-inv_j.getVal(d,0)-inv_j.getVal(d,1)));      
          }          
        }//TRIANGLE
			} else { //!!!DIM 3
          if (m_nodxelem==8){
            for (int d=0;d<m_dim;d++){ //HEXA
              jacob.Set(0,d,0.125*(-x2.getVal(0,d)+x2.getVal(1,d)+x2.getVal(2,d)-x2.getVal(3,d)-x2.getVal(4,d)+x2.getVal(5,d)+x2.getVal(6,d)-x2.getVal(7,d)));  
              jacob.Set(1,d,0.125*(-x2.getVal(0,d)-x2.getVal(1,d)+x2.getVal(2,d)+x2.getVal(3,d)-x2.getVal(4,d)-x2.getVal(5,d)+x2.getVal(6,d)+x2.getVal(7,d)));  
              jacob.Set(2,d,0.125*(-x2.getVal(0,d)-x2.getVal(1,d)-x2.getVal(2,d)-x2.getVal(3,d)+x2.getVal(4,d)+x2.getVal(5,d)+x2.getVal(6,d)+x2.getVal(7,d))); 
              //jacob.Set(0,d,-x2.getVal(0,d) + x2.getVal(1,d) + x2.getVal(2,d) - x2.getVal(3,d));  

            }


          AdjMat(jacob, &inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
          //printf("ADJ J ptr\n");
          //inv_j.Print();          //printf("jacob\n");jacob.Print();
                  
          // jacob.Print();
          ////printf("INV J2 not ptr\n");
          //inv.Print();
          
          //inv.Print();
          
          for (int d=0;d<m_dim;d++){            
            dHxy_detJ_loc.Set(d,0,0.125*(-inv_j.getVal(d,0)-inv_j.getVal(d,1)-inv_j.getVal(d,2)));         
            dHxy_detJ_loc.Set(d,1,0.125*( inv_j.getVal(d,0)-inv_j.getVal(d,1)-inv_j.getVal(d,2)));  
            dHxy_detJ_loc.Set(d,2,0.125*( inv_j.getVal(d,0)+inv_j.getVal(d,1)-inv_j.getVal(d,2)));  
            dHxy_detJ_loc.Set(d,3,0.125*(-inv_j.getVal(d,0)+inv_j.getVal(d,1)-inv_j.getVal(d,2)));             
            dHxy_detJ_loc.Set(d,4,0.125*(-inv_j.getVal(d,0)-inv_j.getVal(d,1)+inv_j.getVal(d,2))); 
            dHxy_detJ_loc.Set(d,5,0.125*( inv_j.getVal(d,0)-inv_j.getVal(d,1)+inv_j.getVal(d,2)));
            dHxy_detJ_loc.Set(d,6,0.125*( inv_j.getVal(d,0)+inv_j.getVal(d,1)+inv_j.getVal(d,2)));
            dHxy_detJ_loc.Set(d,7,0.125*(-inv_j.getVal(d,0)+inv_j.getVal(d,1)+inv_j.getVal(d,2)));
          }
          //dHxy_detJ_loc.Mul(0.125); /////.DO NOT USE THIS!! --- ERRORS ---

          // // elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.125d0    
          } else if (m_nodxelem==4){ //TETRA
            //printf("Element %d\n",e);
            //1 - r - s - t, N2 = r, N3 = s, N4 = t, 
            //dHdrs [h1'r, h2,r]
            //      [h1's,
            //dHdrs [-1,1,0,0]x  X1 Y1 Z1
            //       -1,0,1,0, x X2 Y2 Z2
            //       -1,0,0,1] x X3 Y3 Z3
            //                   x4 Y4 Z4
            //J(0,d) =d1-d4, J(0,1)= 
            //J =dr/dx
            for (int d=0;d<m_dim;d++){
              jacob.Set(0,d,x2.getVal(1,d)-x2.getVal(0,d) ); //d1-d4
              jacob.Set(1,d,x2.getVal(2,d)-x2.getVal(0,d) );            
              jacob.Set(2,d,x2.getVal(3,d)-x2.getVal(0,d) );      
            }
            //USE ADJ TO NOT DIVIDE BY DET
            //dHdr = dH/dr dr/dx
            AdjMat(jacob, &inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
            //printf(" J ptr\n");
            //jacob.Print();
            //printf("ADJ J ptr\n");
            //inv_j.Print();          //printf("jacob\n");jacob.Print();
            //invj ((d,X) x dHdrs [-1,1,0,0;  
            //                     -1,0,1,0;
            //                     -1,0,0,1]
            for (int d=0;d<m_dim;d++){    
              /////ROWS OF INVJ
              dHxy_detJ_loc.Set(d,0,-inv_j.getVal(d,0)-inv_j.getVal(d,1)-inv_j.getVal(d,2));   
              dHxy_detJ_loc.Set(d,1, inv_j.getVal(d,0) );     
              dHxy_detJ_loc.Set(d,2, inv_j.getVal(d,1) );     
              dHxy_detJ_loc.Set(d,3, inv_j.getVal(d,2) );     
            }

          }//TETRA

          
      } // end if  !!!!DIM
      
      m_detJ[offset] = jacob.calcDet();
      //printf("det J %f allocated, offset %d\n",m_detJ[offset],offset);
      //if (m_detJ[offset]>1.0 || m_detJ[offset]<1.0e-3){
        //printf("--------------------- WARNGIN JACOBIAN\n");
        
        //if (e<10){
        ///printf("ELNODES %d %d %d %d \n",m_elnod[nind],m_elnod[nind+1],m_elnod[nind+2],m_elnod[nind+3]);
        //printf("det J %f allocated, offset %d\n",m_detJ[offset],offset);
        //}
        
      //}
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
          
          // dHrs.Set(0,0,-1.0*(1-gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs.Set(1,0,-1.0*(1+gpc[gp][0])*(1.0-gpc[gp][2])); dHrs.Set(2,0,-1.0*(1+gpc[gp][0])*(1.0-gpc[gp][1])); //dh1/d(r,s,t)
          // dHrs.Set(0,1,     (1-gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs.Set(1,1,-1.0*(1+gpc[gp][0])*(1.0-gpc[gp][2])); dHrs.Set(2,1,-1.0*(1-gpc[gp][0])*(1.0-gpc[gp][1])); //dh2/d(r,s,t)
					
					// dHrs.Set(0,2,     (1+gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs.Set(1,2,     (1+gpc[gp][0])*(1.0-gpc[gp][2])); dHrs.Set(2,2,-1.0*(1+gpc[gp][0])*(1.0+gpc[gp][1]));
					// dHrs.Set(0,3,-1.0*(1+gpc[gp][1])*(1.0-gpc[gp][2]));  dHrs.Set(1,3,     (1-gpc[gp][0])*(1.0-gpc[gp][2])); dHrs.Set(2,3,-1.0*(1+gpc[gp][0])*(1.0+gpc[gp][1]));
          
          // dHrs.Set(0,4,-1.0*(1-gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs.Set(1,4,-1.0*(1-gpc[gp][0])*(1.0+gpc[gp][2])); dHrs.Set(2,4,     (1-gpc[gp][0])*(1.0-gpc[gp][1]));
          // dHrs.Set(0,5,     (1-gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs.Set(1,5,-1.0*(1+gpc[gp][0])*(1.0+gpc[gp][2])); dHrs.Set(2,5,     (1+gpc[gp][0])*(1.0-gpc[gp][1]));
          
          // dHrs.Set(0,6,     (1+gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs.Set(1,6,     (1+gpc[gp][0])*(1.0+gpc[gp][2])); dHrs.Set(2,6,     (1+gpc[gp][0])*(1.0+gpc[gp][1]));
          // dHrs.Set(0,7,-1.0*(1+gpc[gp][1])*(1.0+gpc[gp][2]));  dHrs.Set(1,7,     (1-gpc[gp][0])*(1.0+gpc[gp][2])); dHrs.Set(2,7,     (1-gpc[gp][0])*(1.0+gpc[gp][1]));
					          

					// //*jacob = 0.125 * MatMul(*dHrs,*x2);
          // MatMul(*dHrs,*x2,jacob);
          // //printf("x2\n");
          // //m_jacob[e].Print();

          // x2.Print();
          // jacob.Mul(0.125);
          //printf("jacob\n");jacob.Print();
          //jacob.Print();
          // ////printf("Jacobian: \n");jacob.Print();
          // //printf("dHrs\n"); dHrs.Print();
           
          // InvMat(*jacob, inv_j);
          // //printf("inv j\n");inv_j.Print();
          // MatMul(*inv_j,*dHrs,dHxy_detJ_loc);
          
          // //printf("Derivative matrix\n");
          // dHxy_detJ_loc.Print();
         
          
          // m_detJ[offset + m_gp_count * gp] = jacob.calcDet();
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
						// dHrs.Set(0,0,-1.0*(1-gpc[gp][1])); dHrs.Set(0,1,     (1-gpc[gp][1])); dHrs.Set(0,2, 1+gpc[gp][1]);   dHrs.Set(0,3,-1.0*(1+gpc[gp][1]));
						// dHrs.Set(1,0,-1.0*(1-gpc[gp][0])); dHrs.Set(1,1,-1.0*(1+gpc[gp][0])); dHrs.Set(1,2,(1+gpc[gp][0]));  dHrs.Set(1,3,     (1-gpc[gp][0]));
					// }
					// //*jacob = 0.125 * MatMul(*dHrs,*x2);
          // MatMul(*dHrs,*x2,jacob);
          // //printf("jacob\n");
          // //m_jacob[e].Print();

          // x2.Print();
          // jacob.Print();
          // jacob.Mul(0.125);
          
          // //m_detJ[offset + m_gp_count * gp] = det(*jacob);
          
          // // elem%dHrs(e,gp,:,:) =  dHrs(:,:)         
          // // !dHrs(2,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))]         
          // // !dHrs(3,:)=[(1+r(i)), (1-r(i)),-(1-r(i)),-(1+r(i))] 
          // // !print *, "dhrs", dHrs 
          // // !print *, "x2", x2 
          // // elem%jacob(e,gp,:,:) = 0.25*matmul(dHrs,x2)
					// //*jacob = 0.25 * MatMul(*dHrs,*x2);
					
					
					// //jacob.Print();
        
        
        
      // }//dim 2 (gp>1)
    }// end if !!gp ==1

    ///// ALLOCATION
    for (int gp=0;gp<m_gp_count;gp++){
      //Domain_d::setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &v)
      //setDerivative(e,gp,dHxy_detJ_loc
      for (int j=0;j<m_nodxelem;j++){
        int offset = e*(m_nodxelem * m_gp_count) + gp * m_nodxelem;
        ////printf ("Offset %d \n", offset);
        
          //m_dH_detJ_dx[offset + j                 ] = dHxy_detJ_loc.operator()(0,j);
          // m_dH_detJ_dx[offset + j] = dHxy_detJ_loc.getVal(0,j);
          // m_dH_detJ_dy[offset + j] = dHxy_detJ_loc.getVal(1,j); 
          // m_dH_detJ_dz[offset + j] = dHxy_detJ_loc.getVal(2,j);      
          setDerivative(e,gp,0,j,dHxy_detJ_loc.getVal(0,j));
          setDerivative(e,gp,1,j,dHxy_detJ_loc.getVal(1,j));
          if (m_dim ==3)
            setDerivative(e,gp,2,j,dHxy_detJ_loc.getVal(2,j));
          //printf("set der: z n %d %f\n",j, dHxy_detJ_loc.getVal(2,j));
          
      }
    }
          
    //printf("jacob\n");
    //jacob.Print();
    //printf("dHdx x detJ\n");
    //dHxy_detJ_loc.Print();
		//printf("END.\n");
    
    //x2.Free();    inv_j.Free();    jacob.Free();    dHxy_detJ_loc.Free();

    //delete x2; //DEFINED ON EACH BLOCK!
    //  delete inv_j, jacob,dHxy_detJ_loc;
  } // e < elem_colunt
  

}

// dev_t void Domain_d::calcElemJAndDerivatives_Tet_SSRI() {
    // // Coordenadas naturales de los puntos de Gauss
    // const double a = 0.58541020;
    // const double b = 0.13819660;
    // double gpc[4][3] = { {a,b,b}, {b,a,b}, {b,b,a}, {b,b,b} };
    
    // par_loop(e, m_elem_count) {
        // Matrix x2(4,3); // Coordenadas nodales
        
        // // 1. Obtener coordenadas nodales
        // for(int i=0; i<4; i++) {
            // int nid = m_elnod[e*4+i];
            // x2.Set(i,0, x[nid].x);
            // x2.Set(i,1, x[nid].y);
            // x2.Set(i,2, x[nid].z);
        // }
        
        // // 2. Calcular en 4 puntos desviadores
        // for(int gp=0; gp<4; gp++) {
            // Matrix dHrs(3,4); // Derivadas naturales
            
            // // Derivadas de funciones de forma en punto gp
            // double r = gpc[gp][0], s = gpc[gp][1], t = gpc[gp][2];
            // dHrs.Set(0,0, -1.0); dHrs.Set(0,1, 1.0); dHrs.Set(0,2, 0.0); dHrs.Set(0,3, 0.0);
            // dHrs.Set(1,0, -1.0); dHrs.Set(1,1, 0.0); dHrs.Set(1,2, 1.0); dHrs.Set(1,3, 0.0);
            // dHrs.Set(2,0, -1.0); dHrs.Set(2,1, 0.0); dHrs.Set(2,2, 0.0); dHrs.Set(2,3, 1.0);
            
            // // Jacobiano y derivadas cartesianas
            // Matrix jacob = MatMul(dHrs, x2);
            // Matrix inv_j = Invert(jacob);
            // Matrix dHxy = MatMul(inv_j, dHrs);
            
            // // Almacenar derivadas para este punto GP
            // for(int n=0; n<4; n++) {
                // setDerivative(e,gp,0,n, dHxy(0,n));
                // setDerivative(e,gp,1,n, dHxy(1,n));
                // setDerivative(e,gp,2,n, dHxy(2,n));
            // }
            // m_detJ[e*5 + gp] = Det(jacob); // 4 primeros puntos
        // }
        
        // // 3. Punto volumétrico (centroide)
        // Matrix dHrs_vol(3,4);
        // dHrs_vol.SetAll(-0.25); // Valor constante en centroide
        
        // Matrix jacob_vol = MatMul(dHrs_vol, x2);
        // m_detJ[e*5 + 4] = Det(jacob_vol); // 5to punto
    // }
// }

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
  return 1;
}
  
dev_t void Domain_d::Calc_Element_Radius() //For axisymm
{
 /*
  do e=1, elem_count
    do i=1,nodxelem
        !print *, "elnod " , elem%elnod(e,i)
        x2(i,:)=nod%x(elem%elnod(e,i),:)
    end do
    do gp=1, elem%gausspc(e)
      !if (bind_dom_type .eq. 3) then 
        elem%radius(e,gp)= DOT_PRODUCT (elem%math(e,gp, 1,:), x2(:,1))
        !print *, "radius", elem%radius(e,gp)
        ! if (axisymm_vol_weight) then
          ! elem%detJ(e,:) = elem%detJ(e,:) * radius
        ! end if
      !end if 
    end do 
  end do
  */
  //// NEED TO KNOW MASS MATRICES
  for (int e=0;e<m_elem_count;e++){
    //double temp = 0.0;
    m_radius[e] = 0.0;
    int offset = m_nodxelem * e;
    for (int ne=0; ne<m_nodxelem;ne++){
       int n = elnod_h[offset+ne];
    //for (int gp=0;gp<m_gp_count;gp++){    
      
      //!if (bind_dom_type .eq. 3) then 
      //  elem%radius(e,gp)= DOT_PRODUCT (elem%math(e,gp, 1,:), x2(:,1))    
      //TO MODIFY BY H matrix
      if (m_dim ==2)
        m_radius[e] += getPosVec2(n).x;
      else
         m_radius[e] += getPosVec3(n).x;
    }
    
    for (int gp=0;gp<m_gp_count;gp++)
      m_radius[m_gp_count*e+gp] /= m_nodxelem;
  } 
  
  
}

dev_t void Domain_d::printVec(double *v){
  //printf("Printing vectors: \n");
  for (int n=0;n<m_node_count;n++){
    for (int d=0;d<m_dim;d++)
      //printf("%.6e ",v[m_dim*n + d]);
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


//// TODO: SOLVE THIS EITHER WITH DOUBLE2 include
#ifndef CUDA_BUILD
inline double length(const double2 &v) {
    return sqrt(v.x * v.x + v.y * v.y);
}
#endif

////// NEW INCLUDING ANGLES
dev_t void Domain_d::calcMinEdgeLength() {
    double min_len = 1.0e6;
    double min_height = 1.0e6;
    double min_angle = 1.0e6;
    double max_angle = -1.0e6;
    m_min_Jnorm = 1.0e6;
    double Jnorm;
    m_bad_elem_count = 0;
    for (int e = 0; e < m_elem_count; e++) {
        double elem_min_height = 1.0e6;
        double elem_min_angle = 1.0e6;
        double elem_max_angle = -1.0e6;

        int off = m_nodxelem * e;

        if (m_dim == 3) {  // --- Tetrahedron case ---
            int a = m_elnod[off];
            int b = m_elnod[off + 1];
            int c = m_elnod[off + 2];
            int d = m_elnod[off + 3];
          
            double3 A = getPosVec3(a);
            double3 B = getPosVec3(b);
            double3 C = getPosVec3(c);
            double3 D = getPosVec3(d);

            // --- edge lengths (igual que antes) ---
            double3 edges[6] = {B-A, C-A, D-A, C-B, D-B, D-C};
            for (int i = 0; i < 6; i++) {
                double len = length(edges[i]);
                if (len < min_len) min_len = len;
            }

            // --- heights (igual que antes) ---
            double3 faces[4][3] = {{B, C, D}, {A, C, D}, {A, B, D}, {A, B, C}};
            for (int i = 0; i < 4; i++) {
                double3 x1 = faces[i][1] - faces[i][0];
                double3 x2 = faces[i][2] - faces[i][0];
                double3 normal = cross(x1, x2);
                double area = length(normal);
                if (area < 1e-12) continue;
                normal = normal / area;

                double3 vec = getPosVec3(m_elnod[off + i]) - faces[i][0];
                double height = fabs(dot(vec, normal));
                if (height < elem_min_height) elem_min_height = height;
            }
            m_elem_length[e] = elem_min_height;
            if (elem_min_height < min_height) min_height = elem_min_height;

            // --- angle computation (dihedral angles between faces) ---
            int face_nodes[4][3] = {
                {b, c, d}, {a, c, d}, {a, b, d}, {a, b, c}
            };
            double3 normals[4];
            for (int i = 0; i < 4; i++) {
                double3 P1 = getPosVec3(face_nodes[i][0]);
                double3 P2 = getPosVec3(face_nodes[i][1]);
                double3 P3 = getPosVec3(face_nodes[i][2]);
                double3 n = cross(P2 - P1, P3 - P1);
                double ln = length(n);
                if (ln > 1e-16) n = n / ln;
                normals[i] = n;
            }
            // ANGLES OLD WAY
            // 6 dihedral angles between face pairs
            // int pairs[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
            // for (int p = 0; p < 6; p++) {
                // double cang = dot(normals[pairs[p][0]], normals[pairs[p][1]]);
                // cang = fmax(-1.0, fmin(1.0, cang));
                // double ang_deg = acos(cang) * 180.0 / M_PI;
                // if (ang_deg < elem_min_angle) elem_min_angle = ang_deg;
                // if (ang_deg > elem_max_angle) elem_max_angle = ang_deg;
            // }
            
            
            
            ////ANGLES; NEW WAY
            //////////////////////--- dihedral angles per edge (robust) ---
            // int edges[6][2] = {
                // {a,b}, {a,c}, {a,d},
                // {b,c}, {b,d}, {c,d}
            // };

            // For each edge, the two opposite faces
            // Each face is given by the 3 nodes that define it
            int faces_per_edge[6][2][3] = {
                {{a,b,c},{a,b,d}}, // edge ab
                {{a,c,b},{a,c,d}}, // edge ac
                {{a,d,b},{a,d,c}}, // edge ad
                {{b,c,a},{b,c,d}}, // edge bc
                {{b,d,a},{b,d,c}}, // edge bd
                {{c,d,a},{c,d,b}}  // edge cd
            };

            for (int e2 = 0; e2 < 6; e2++) {

                // --- face 1 normal ---
                double3 P1 = getPosVec3(faces_per_edge[e2][0][0]);
                double3 P2 = getPosVec3(faces_per_edge[e2][0][1]);
                double3 P3 = getPosVec3(faces_per_edge[e2][0][2]);
                double3 n1 = cross(P2 - P1, P3 - P1);
                double l1 = length(n1);
                if (l1 < 1e-16) continue;
                n1 = n1/l1;

                // --- face 2 normal ---
                P1 = getPosVec3(faces_per_edge[e2][1][0]);
                P2 = getPosVec3(faces_per_edge[e2][1][1]);
                P3 = getPosVec3(faces_per_edge[e2][1][2]);
                double3 n2 = cross(P2 - P1, P3 - P1);
                double l2 = length(n2);
                if (l2 < 1e-16) continue;
                n2 = n2/l2;

                // --- interior dihedral angle ---
                double cang = dot(n1, n2);
                cang = fmax(-1.0, fmin(1.0, cang));
                double ang = M_PI - acos(cang);   // interior angle

                double ang_deg = ang * 180.0 / M_PI;

                elem_min_angle = fmin(elem_min_angle, ang_deg);
                elem_max_angle = fmax(elem_max_angle, ang_deg);
            }
            
            double3 BA = B - A;
            double3 CA = C - A;
            double3 DA = D - A;
            
            // --- edge lengths ---
            double l2_sum = 0.0;
            double3 E[6] = {B-A, C-A, D-A, C-B, D-B, D-C};
            for (int i = 0; i < 6; i++) {
                double l = length(E[i]);
                l2_sum += l*l;
            }
            //isthe same that m_detJ but they are calculated after
            double Jgeom = dot(BA, cross(CA, DA));  // = 6*V (signed)
            //cout << "Jgeom "<<Jgeom<<", detJ "<<m_detJ[e]<<endl;

            double V = fabs(Jgeom) / 6.0;
            Jnorm = (6.0 * V) / pow(l2_sum, 1.5);
            
            if (Jnorm<m_min_Jnorm) m_min_Jnorm = Jnorm;
            

        } else if (m_dim == 2) {  // --- Triangle case ---
            int a = m_elnod[off];
            int b = m_elnod[off + 1];
            int c = m_elnod[off + 2];

            double2 A = getPosVec2(a);
            double2 B = getPosVec2(b);
            double2 C = getPosVec2(c);

            // --- edges ---
            double2 AB = B - A;
            double2 BC = C - B;
            double2 CA = A - C;

            double lenAB = length(AB);
            double lenBC = length(BC);
            double lenCA = length(CA);
            if (lenAB < min_len) min_len = lenAB;
            if (lenBC < min_len) min_len = lenBC;
            if (lenCA < min_len) min_len = lenCA;

            // --- area/height (igual que antes) ---
            double area = 0.5 * fabs((B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y));
            if (area > 1e-12) {
                double base = lenAB;
                double height = 2.0 * area / base;
                elem_min_height = height;
                if (height < min_height) min_height = height;
            }
            m_elem_length[e] = elem_min_height;

            // --- internal angles (ley de cosenos) ---
            double cosA = (lenAB*lenAB + lenCA*lenCA - lenBC*lenBC) / (2.0 * lenAB * lenCA);
            double cosB = (lenAB*lenAB + lenBC*lenBC - lenCA*lenCA) / (2.0 * lenAB * lenBC);
            double cosC = (lenBC*lenBC + lenCA*lenCA - lenAB*lenAB) / (2.0 * lenBC * lenCA);
            cosA = fmax(-1.0, fmin(1.0, cosA));
            cosB = fmax(-1.0, fmin(1.0, cosB));
            cosC = fmax(-1.0, fmin(1.0, cosC));
            double angA = acos(cosA) * 180.0 / M_PI;
            double angB = acos(cosB) * 180.0 / M_PI;
            double angC = acos(cosC) * 180.0 / M_PI;
            elem_min_angle = fmin(angA, fmin(angB, angC));
            elem_max_angle = fmax(angA, fmax(angB, angC));
        }

        //cout << "angle element "<<e<<endl;
        //m_elem_min_angle[e] = elem_min_angle;
        //m_elem_max_angle[e] = elem_max_angle;
        
        if (elem_min_angle < min_angle) min_angle = elem_min_angle;
        if (elem_max_angle > max_angle) max_angle = elem_max_angle;
    
    
      if (elem_min_angle<5.0 || Jnorm<0.001)
        m_bad_elem_count++;
      
    }//elem e
    
    // cout <<"elem_min_angle"<<min_angle<<endl;
    // cout <<"elem_max_angle"<<max_angle<<endl;
    // cout <<"Min Jnorm "<<m_min_Jnorm<<endl;
    
    // global min/max
    m_min_length = min_len;
    m_min_height = min_height;
    m_min_angle  = min_angle;
    m_max_angle  = max_angle;
}


///// ALREADY ALLOCATED
void Domain_d::setNode(const int &i, const double &_x, const double &_y, const double &_z){
  if (i<m_node_count){
  x[i*m_dim  ]=_x;
  x[i*m_dim+1]=_y;
  x[i*m_dim+2]=_z;
  //return 1;
  }
else{cout << "Node allocation error, node pos larger than node count."<<endl;}
      //return 0;
}
  
dev_t void Domain_d::BlendStresses(const double &s, const double &pl_strain_max){
    

      //double s = pow((step_count - last_step_remesh)/(double)STEP_RECOV, 2.0); // Curva cuadrática
      
      // 1. Stress Blending con conservación de energía
      for (int e=0; e<m_elem_count; ++e) {
          //double weight = s * std::min(1.0, pl_strain[e]/pl_strain_max); // Peso adaptativo
          double weight = s;
          for (int c=0; c<6; ++c) {
              m_sigma[e*6 + c] = 
                  weight * m_sigma[e*6 + c] + 
                  (1.0 - weight) * m_sigma_prev[e*6 + c];

              //~ m_tau[e*6 + c] = 
                  //~ weight * m_sigma[e*6 + c] + 
                  //~ (1.0 - weight) * m_tau_prev[e*6 + c];
                  
          }
      }
      
      // 2. Plastic Strain Blending no-lineal (evita discontinuidades)
      for (int e=0; e<m_elem_count; ++e) {
          //double alpha = s * (1.0 - exp(-5.0*pl_strain[e]/pl_strain_ref));
          double alpha = s;
          //double alpha = s * (1.0 - exp(-5.0 * pl_strain[e] / pl_strain_max));
          pl_strain[e] = alpha * pl_strain[e] + (1.0 - alpha) * pl_strain_prev[e];
          
          // Asegurar monotonía para modelos plásticos
          if (pl_strain[e] < pl_strain_prev[e]) {
              pl_strain[e] = pl_strain_prev[e] + 1e-6*(pl_strain[e] - pl_strain_prev[e]);
          }
      }


}

dev_t void Domain_d::BlendField(const double &s, const int size, const int &d, double *prev, double *curr) {
    

      //double s = pow((step_count - last_step_remesh)/(double)STEP_RECOV, 2.0); // Curva cuadrática
      
      // 1. Stress Blending con conservación de energía
      for (int e=0; e<size; ++e) {
          //double weight = s * std::min(1.0, pl_strain[e]/pl_strain_max); // Peso adaptativo
          double weight = s;
          for (int c=0; c<d; ++c) {
              curr[e*6 + c] = 
                  weight * prev[e*d + c] + 
                  (1.0 - weight) * prev[e*d + c];
                  
          }
      }

}

dev_t void Domain_d::postRemeshGlobFilter()

    //const double* m_vprev,       // Velocidades PRE-remallado [nnode*m_dim]
    //const double* pl_strain,     // Deformación plástica nodal [nnode]
    //const double* nodal_mass,    // Masa nodal (opcional, para ponderar)
{
    double kv=0.6;
    double ka=0.2;
    double strain_threshold=0.1;

    // ---- 1. Calcular promedio ponderado de velocidad ----
    double v_avg[3] = {0,0,0}, mass_total=0;
    for(int i=0; i<m_node_count; i++) {
        if(pl_strain[i] < strain_threshold) continue;
        //double w = nodal_mass ? nodal_mass[i] : 1.0;
        double w = m_mdiag[i];
        for(int d=0; d<m_dim; d++) 
            v_avg[d] += w * m_vprev[m_dim*i+d]; // Usamos m_vprev como referencia!
        mass_total += w;
    }
    for(int d=0; d<m_dim; d++) 
        v_avg[d] /= (mass_total > 0 ? mass_total : 1.0);

    // ---- 2. Filtrado conservativo ----
    for(int i=0; i<m_node_count; i++) {
        if(pl_strain[i] < strain_threshold) continue;

        // Magnitudes pre/post
        double v_mag_prev = 0.0, v_mag_current = 0.0;
        for(int d=0; d<m_dim; d++) {
            v_mag_prev += m_vprev[m_dim*i+d] * m_vprev[m_dim*i+d];
            v_mag_current += v[m_dim*i+d] * v[m_dim*i+d];
        }
        v_mag_prev = sqrt(v_mag_prev);
        v_mag_current = sqrt(v_mag_current);

        // Factor de blending basado en plasticidad
        #ifndef CUDA_BUILD
        double alpha = std::min(1.0, pl_strain[i]/(2.0*strain_threshold));
        #else
        double alpha = min_(1.0, pl_strain[i]/(2.0*strain_threshold));
        #endif
        double v_target = (1.0-alpha)*v_mag_prev + alpha*v_mag_current;

        // Aplicar filtro (kv=1: no filtro, kv=0: full damping)
        for(int d=0; d<m_dim; d++) {
            // Componente fluctuación respecto al promedio
            double v_fluc = v[m_dim*i+d] - v_avg[d];
            v[m_dim*i+d] = v_avg[d] + kv * v_fluc;

            // Limitador físico basado en m_vprev
            double dv_max = 0.1 * v_mag_prev + 1e-8;
            if(fabs(v[m_dim*i+d] - m_vprev[m_dim*i+d]) > dv_max) {
                v[m_dim*i+d] = m_vprev[m_dim*i+d] + 
                    (v[m_dim*i+d] > m_vprev[m_dim*i+d] ? dv_max : -dv_max);
            }

            // Filtrado de aceleración (más agresivo)
            //a[m_dim*i+d] *= ka;
        }
    }
}

dev_t void Domain_d::SmoothDeviatoricStress(double alpha) {
    // Suavizado Laplaciano de τ (conservativo)
    double* tau_smoothed = new double[m_elem_count * 6];
    #pragma omp parallel for
    for (int e = 0; e < m_elem_count; ++e) {
        tensor3 tau_avg;
        clear(tau_avg);
        int n_neighbors = 0;
        for (int a = 0; a < m_nodxelem; ++a) {
            int nid = m_elnod[e*m_nodxelem + a];
            for (int i = 0; i < m_nodel_count[nid]; ++i) {
                int e_neigh = m_nodel[m_nodel_offset[nid] + i];
                if (e_neigh != e) {
                    tau_avg = tau_avg + FromFlatSym(m_tau, e_neigh * 6);
                    n_neighbors++;
                }
            }
        }
        tau_avg = (n_neighbors > 0) ? tau_avg * (1.0 / n_neighbors) : FromFlatSym(m_tau, e * 6);
        // Blending suave
        tensor3 tau_e = FromFlatSym(m_tau, e * 6);
        tau_e = alpha * tau_avg + (1.0 - alpha) * tau_e;
        ToFlatSymPtr(tau_e, tau_smoothed, e * 6);
    }
    memcpy(m_tau, tau_smoothed, sizeof(double) * m_elem_count * 6);
    delete[] tau_smoothed;
}

dev_t double Domain_d::getPtrMax(double *v, const int &size, const int &dim){
  double max=0.0;
  for (int i=0;i<size;i++){
    double norm2;
    if (dim==3){
      vector_t val;
      val.x=v[m_dim*i];      val.y=v[m_dim*i+1];      val.z=v[m_dim*i+2];
      norm2=val.x*val.x+val.y*val.y+val.z*val.z;
    }
    if (norm2>max)max=norm2;
  }
  return max;
}

void Domain_d::CorrectLocalVelocityPeaks() {
    const double v_ref = 1.2;  // Tu velocidad característica [m/s]
    const double v_limit = 10.0 * v_ref; // Umbral para corrección (ej. 12 m/s)
    int corrected_nodes = 0;

    // --- Paso 1: Identificar nodos problemáticos ---
    #pragma omp parallel for reduction(+:corrected_nodes)
    for (int n = 0; n < m_node_count; n++) {
        vector_t vel = getVelVec(n);
        double v_mag = norm(vel);
        
        if (v_mag > v_limit) {
            // --- Paso 2: Corrección física conservando dirección ---
            double correction_factor = v_limit / v_mag;
            
            for (int d = 0; d < m_dim; d++) {
                v[m_dim * n + d] *= correction_factor;
                a[m_dim * n + d] *= 0.5 * correction_factor; // Amortiguar aceleración también
            }
            corrected_nodes++;
        }
    }
  }

    // dev_t void Domain::Save_Step() {
        // // Copiar tensores elementales
        // memcpy(m_sigma_prev, m_sigma, sizeof(double) * m_elem_count * 6);
        // memcpy(m_str_rate_prev, m_str_rate, sizeof(double) * m_elem_count * 6);
        // memcpy(m_tau_prev, m_tau, sizeof(double) * m_elem_count * 6);
        
        // // Copiar variables nodales
        // memcpy(m_u_prev, u, sizeof(double) * m_node_count * m_dim);
        // memcpy(m_v_prev, v, sizeof(double) * m_node_count * m_dim);
        // memcpy(m_a_prev, a, sizeof(double) * m_node_count * m_dim);
        // memcpy(m_prev_a_prev, prev_a, sizeof(double) * m_node_count * m_dim);
        
        // // Copiar variables escalares elementales
        // memcpy(m_pl_strain_prev, pl_strain, sizeof(double) * m_elem_count);
        // memcpy(m_vol_prev, vol, sizeof(double) * m_elem_count);
        // memcpy(m_rho_prev, rho, sizeof(double) * m_elem_count);
        // memcpy(m_p_prev, p, sizeof(double) * m_elem_count);
    // }
    
    // // Función para restaurar el estado anterior
    // dev_t void Domain_d::Restore_Step() {
        // // Restaurar tensores elementales
        // memcpy(m_sigma, m_sigma_prev, sizeof(double) * m_elem_count * 6);
        // memcpy(m_str_rate, m_str_rate_prev, sizeof(double) * m_elem_count * 6);
        // memcpy(m_tau, m_tau_prev, sizeof(double) * m_elem_count * 6);
        
        // // Restaurar variables nodales
        // memcpy(u, m_u_prev, sizeof(double) * m_node_count * m_dim);
        // memcpy(v, m_v_prev, sizeof(double) * m_node_count * m_dim);
        // memcpy(a, m_a_prev, sizeof(double) * m_node_count * m_dim);
        // memcpy(prev_a, m_prev_a_prev, sizeof(double) * m_node_count * m_dim);
        
        // // Restaurar variables escalares elementales
        // memcpy(pl_strain, m_pl_strain_prev, sizeof(double) * m_elem_count);
        // memcpy(vol, m_vol_prev, sizeof(double) * m_elem_count);
        // memcpy(rho, m_rho_prev, sizeof(double) * m_elem_count);
        // memcpy(p, m_p_prev, sizeof(double) * m_elem_count * m_gp_count);
    // }


//~ void post_remesh_accel_filter(
    //~ double* a, 
    //~ int nnode, 
    //~ int dim, 
    //~ double ka=0.2,  // Factor de filtrado (ajustable)
    //~ bool is_first_step=false
//~ ) {
    //~ if (is_first_step) ka = 0.2;  // Filtrado más agresivo al inicio

    //~ // Filtro pasa-bajas simple (promedio con vecinos)
    //~ #pragma omp parallel for
    //~ for (int i=0; i<nnode; ++i) {
        //~ double a_avg[3] = {0};
        //~ int count = 0;
        //~ for (int j : neighbors[i]) {  // Asume vecindario definido
            //~ for (int d=0; d<dim; ++d) a_avg[d] += a[j*dim + d];
            //~ count++;
        //~ }
        //~ for (int d=0; d<dim; ++d) {
            //~ a[i*dim + d] = ka * a[i*dim + d] + (1.0 - ka) * (a_avg[d] / (count + 1e-16));
        //~ }
    //~ }
//~ }

double3 Domain_d::ComputeCentroid() const {
    // Verificación básica
    if (m_node_count == 0 || x == nullptr) {
        std::cerr << "Error: Domain_d::ComputeCentroid called on empty domain." << std::endl;
        return make_double3(0.0, 0.0, 0.0);
    }

    double3 centroid = make_double3(0.0, 0.0, 0.0);

    // Sumar todas las posiciones nodales
    for (int n = 0; n < m_node_count; ++n) {
        centroid.x += x[3 * n + 0];
        centroid.y += x[3 * n + 1];
        centroid.z += x[3 * n + 2];
    }

    // Promediar
    centroid.x /= m_node_count;
    centroid.y /= m_node_count;
    centroid.z /= m_node_count;

    return centroid;
}


void Domain_d::setFixSymm(){
  double symtol = 1.0e-4;

  for (int i=0;i<getNodeCount();i++){
    for (int d=0;d<3;d++){
      if (m_symm[d]){
        #ifdef CUDA_BUILD
        #else
          double coord;
          if      (d==0)  coord = getPosVec3(i).x;
          else if (d==1)  coord = getPosVec3(i).y;
          else            coord = getPosVec3(i).z;
          //cout << "BCZ COUNT "<<bcz_nod_h.size()<<endl;
          if (coord < symtol ) {
            AddBCVelNode(i,d,0);
          }
        #endif
      }
    }
  }
}

void Domain_d::addMeshData(const TriMesh_d &m){
  trimesh->AddMesh(m);
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

__global__ void calcMinEdgeLength(Domain_d *dom_d){
  dom_d->calcMinEdgeLength();
}

__global__ void calcMinEdgeLengthKernel(Domain_d *dom_d){
  dom_d->calcMinEdgeLength();
}

__global__ void InitElemValuesKernel(Domain_d *dom_d, double *arr, double val){
  dom_d->InitElemValues(arr, val);
}

__global__ void InitStressesFromMatKernel(Domain_d *dom_d){
  dom_d->InitStressesFromMat();
 }

  // Kernel to retrieve only the private value
__global__ void getMinLengthKernel(Domain_d *dom_d, double *d_value) {
    *d_value = dom_d->getMinLength();  // Store private value in separate memory
}

};
	
