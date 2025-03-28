#ifndef _REMESHER_
#define _REMESHER_

namespace MetFEM{

class Domain_d;

class ReMesher{
  public:
  ReMesher(){}
  ReMesher(Domain_d *);
  
  void Map();
  
  void ReMesh();

  protected:
  Domain_d *m_dom;
};

}; //MetFEM
#endif
