#ifndef _REMESHER_
#define _REMESHER_

#ifdef REMESH_OMEGA_H
#include <Omega_h_mesh.hpp>
#endif

#include <vector>
#include <array>
#include <cmath>


struct Mesh{
  
  
};


enum Remesh_Type { MMG=0, OMG_H=1 };

#ifdef REMESH_OMEGA_H
void create_mesh(Omega_h::Mesh& mesh, 
#ifdef CUDA_BUILD
                 double* d_node_coords, int num_nodes, 
                 int* d_element_conn, int num_elements
#else
                 double* h_node_coords, int num_nodes, 
                 int* h_element_conn, int num_elements
#endif
                 ) ;
                 
#else
void create_mesh(Mesh& mesh, 
                 double* h_node_coords, int num_nodes, 
                 int* h_element_conn, int num_elements  
);
#endif
                 
                 
#include "defs.h"

namespace MetFEM{

class Domain_d;

struct remesh_params {
    double nodal_tol;   //fract 
    double lambda_tol;  //element inside
    bool   map_momentum;
    
}; 

class ReMesher{
  public:
  
  Remesh_Type m_type;
  
  ReMesher()//:dim(3)
  {}
  ReMesher(Domain_d *);
  
  void Generate_omegah();
  void Generate_mmg();
  
  //dest,origin,switches between mesher
  void MapNodal(double *vfield, double *o_field); //VECTORIAL @TODO: ADD TEMPLATE
  void MapNodalScalar(double *vfield, double *o_field);
  void MapElem(double *vfield, double *o_field, int field_dim=1);
  
  //THIS USES m_x,m_elnod,m_xx_count as output mesh data
  template <int dim>
  void MapNodalVectorRaw(double *vfield, double *o_field); //
  void MapElemVectorRaw  (double *vfield, double *o_field, int field_dim); ///

  void MapElemVectorL2(double *vfield, double *o_field, int field_dim);
  
  ///// THIS IS FOR MIN SQUARES EIGEN BASED ////////////

  //~ std::vector<int> ReMesher::find_nearest_old_elems(const std::array<double, 3>& point, 
                                                   //~ int num_neighbors,
                                                   //~ const std::vector<std::array<double, 3>>& old_centroids);
  //~ std::array<double, 3> ReMesher::get_old_elem_centroid(int old_elem);                                                   
  //~ std::array<double, 3> ReMesher::get_new_elem_centroid(int new_elem) ;
  
  //////////////////////////////////// NO EIGEN BASED VERSION
    std::array<double, 3> get_old_elem_centroid(int old_elem);
    std::array<double, 3> get_new_elem_centroid(int new_elem);
    std::vector<int> find_nearest_old_elems(const std::array<double, 3>& point, 
                                          int num_neighbors,
                                          const std::vector<std::array<double, 3>>& old_centroids);
    bool solve_4x4_system(double A[4][4], double b[4], double x[4]);
    double distance(const std::array<double, 3>& a, const std::array<double, 3>& b);
    
    // Función de búsqueda de K vecinos más cercanos
    std::vector<int> find_k_closest_nodes(const std::array<double, 3>& point, int k);
    int find_closest_node(const std::array<double, 3>& point); // Sobrecarga para array
  //////////////////////////////////////////
  
  void FindMapElemClosest();
  
  #ifdef REMESH_OMEGA_H
  template <int dim>
  void MapNodalVector(Omega_h::Mesh &mesh, double *, double *);
  template <int dim>
  void MapElemVector (Omega_h::Mesh &mesh, double *, double *, int field_dim=1);
  template <int dim>
  void ProjectElemToNodes(Omega_h::Mesh& mesh,  double* elemvals, double* nodal_vals, int field_dim);
  template <int dim>
  void ProjectNodesToElem(Omega_h::Mesh& mesh,  double* nodal_vals, double* elemvals, int field_dim);
  template <int dim>
  void HybridProjectionElemToElem(Omega_h::Mesh& mesh, double* new_elemvals,  double* old_elemvals, int field_dim=1);
  template <int dim, typename T>
  void MapElemPtrVector(Omega_h::Mesh& mesh, T* vfield, T* o_field);
  #endif
  
  int find_closest_node(const double x[3]);
  
  void ReMapBCsByFace(int* old_bc_nod,
                        double* old_bc_val,
                        int* new_bc_nod,
                        double* new_bc_val,
                        int bc_count);

  template <int dim>
  void MapElemVectors ();
  

  
  void WriteDomain();
      
  void MeshMMG(); 
  void ReMesh();
  void ReMapBCs(int  *old_bc_nod,
                      double *old_bc_val,

                    int  *new_bc_nod,
                    double *new_bc_val,
                    int bc_count);
                    
  ~ReMesher(){
    delete[] m_x,m_elnod;//CONVERT TO free to use in CUDA
    free_t(m_closest_elem);
  }
  
  remesh_params m_params;
  
  protected:
  
  

  Domain_d *m_dom;
  int *m_mapelem; //DEVICE VECTOR, MAPS ORIGINAL VECTOR TO NEW MESH VECTOR
  #ifdef REMESH_OMEGA_H 
  Omega_h::Mesh m_mesh;
  Omega_h::Mesh m_old_mesh;
  #else
  Mesh m_mesh;
  Mesh m_old_mesh;
  #endif
  //IF OMEGA H NOT USED
  double *m_x;
  int    *m_elnod;
  int     m_node_count;
  int     m_elem_count;
  
  int *m_closest_elem;
  
  ////// FOR MIN SQ 
  std::vector<std::array<double, 3>> old_node_coords; // Coordenadas de nodos OLD
  
  
};

}; //MetFEM
#endif
