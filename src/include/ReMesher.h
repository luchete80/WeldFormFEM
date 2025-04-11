#ifndef _REMESHER_
#define _REMESHER_

#include <Omega_h_mesh.hpp>

void create_mesh(Omega_h::Mesh& mesh, 
#ifdef CUDA_BUILD
                 double* d_node_coords, int num_nodes, 
                 int* d_element_conn, int num_elements
#else
                 double* h_node_coords, int num_nodes, 
                 int* h_element_conn, int num_elements
#endif
                 ) ;

namespace MetFEM{

class Domain_d;

class ReMesher{
  public:
  ReMesher()//:dim(3)
  {}
  ReMesher(Domain_d *);
  
  void Generate_omegah();
  void Generate_mmg();
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


  template <int dim>
  void MapElemVectors ();
  
  template <int dim, typename T>
  void MapElemPtrVector(Omega_h::Mesh& mesh, T* vfield, T* o_field);
  
  void WriteDomain();
      
  void MeshMMG(); 
  void ReMesh();
  void ReMapBCs(int  *old_bc_nod,
                      double *old_bc_val,

                    int  *new_bc_nod,
                    double *new_bc_val,
                    int bc_count);
  protected:

  Domain_d *m_dom;
  int *m_mapelem; //DEVICE VECTOR, MAPS ORIGINAL VECTOR TO NEW MESH VECTOR
  Omega_h::Mesh m_mesh;
  Omega_h::Mesh m_old_mesh;
};

}; //MetFEM
#endif
