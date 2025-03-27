#ifndef _REMESHER_
#define _REMESHER_

#include <Omega_h_mesh.hpp>


namespace MetFEM{

class Domain_d;

class ReMesher{
  public:
  ReMesher(){}
  ReMesher(Domain_d *);
  

  template <int dim>
  void Map(Omega_h::Mesh &mesh);
  
  void MeshMMG(); 
  void ReMesh();

  protected:
  Domain_d *m_dom;
  Omega_h::Mesh m_mesh;
};

}; //MetFEM
#endif
