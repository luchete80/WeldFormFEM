
#include "Domain_d.h"
#include <string>
#include <iostream>

using namespace std;
using namespace LS_Dyna;

namespace MetFEM{

Domain_d::Domain_d(std::string name){
  cout << "Reading "<<name<<endl;
  lsdynaReader reader(name.c_str());
  
  vector_t *x_H =  new vector_t [reader.m_node.size()];
  
  for (int n=0;n<reader.m_node.size();n++){
    //reader.m_node[n].m_x[0]
  }

}

};
