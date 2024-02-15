#include "lsdynaReader.h"
#include <iostream>

using namespace std;

int main(){
  
  lsdynaReader reader("sphere-plate.k");
  cout << "---------------" << endl;
  cout << "Node    count: " << reader.m_node.size() <<endl;
  cout << "Element count: " << reader.m_elem.size() <<endl;
  cout << "BOUNDARY_SPC_NODE count: " << reader.m_spc_nod.size()<<endl;
  
  //TODO: 
  //READ ALL BOUNDARY_SPC_NODE instances
  // *BOUNDARY_SPC_SET
  // $#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
           // 1         0         1         1         1         1         1         1
  // *SET_NODE_LIST_TITLE
  // *BOUNDARY_PRESCRIBED_MOTION_NODE
  //BOUNDARY_PRESCRIBED_MOTION_SET
}