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
  

  template <int dim>
  void MapNodalVector(Omega_h::Mesh &mesh, double *, double *);
  template <int dim>
  void MapElemVector (Omega_h::Mesh &mesh, double *, double *, int field_dim=1);
  template <int dim, typename T>
  void MapElemPtrVector(Omega_h::Mesh& mesh, T* vfield, T* o_field);
  
  void WriteDomain();
      
  void MeshMMG(); 
  void ReMesh();

  protected:

  Domain_d *m_dom;
  Omega_h::Mesh m_mesh;
};

}; //MetFEM
#endif
