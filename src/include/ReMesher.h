#ifndef _REMESHER_
#define _REMESHER_

#include <Omega_h_mesh.hpp>


namespace MetFEM{

class Domain_d;

class ReMesher{
  public:
  ReMesher()//:dim(3)
  {}
  ReMesher(Domain_d *);
  
  void Generate_mmg();
  template <int dim>
  void MapNodalVector(Omega_h::Mesh &mesh, double *, double *);
  template <int dim>
  void MapElemVector (Omega_h::Mesh &mesh, double *, double *, int field_dim=1);

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
