/*************************************************************************/
/*  inpout.C                                                     */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                              */
/*  https://www.opensourcemech.com                                                                */
/*                                                                       */
/*  Copyright (c) 2025-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/

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
