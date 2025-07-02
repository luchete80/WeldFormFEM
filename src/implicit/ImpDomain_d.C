#include "Domain_d.h"
#include <iostream>
#include <vector>

#include "Matrix.h"
#include <sstream>
#include <fstream> 
#include <iostream>


#include "Tensor3.C"
#include "lsdynaReader.h"


using namespace std;
using namespace LS_Dyna;

namespace MetFEM {
  
#define ELNOD  4   //ORIGINALLY 8
#define FACENOD 3  //ORIGINALLY 4
#define ELFAC  4   //ORIGINALLY 6
//OLD FOT HEXA, CHANGE IT





void Domain_d::AssignMaterial (Material_ *material_h) {
//    cudaMalloc((void**)&materials, 1 * sizeof(Material_ )); //
    printf("Assigning Material\n");
    malloc_t(materials, Material_,1);
    memcpy_t(materials, material_h, 1 * sizeof(Material_));	
}


  ///// ALREADY ALLOCATED
  //~ void Domain_d::setNode(const int &i, const double &_x, const double &_y, const double &_z){
    //~ if (i<m_node_count){
    //~ x[i*m_dim  ]=_x;
    //~ x[i*m_dim+1]=_y;
    //~ x[i*m_dim+2]=_z;
    //~ //return 1;
    //~ }
  //~ else{cout << "Node allocation error, node pos larger than node count."<<endl;}
        //~ //return 0;
  //~ }

};
