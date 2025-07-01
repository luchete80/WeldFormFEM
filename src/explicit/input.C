/***********************************************************************************
* WeldFormFEM - A C++/CUDA library to simulate Mechanical Systems using            *
*               explicit Finite Element Methods                                    *
* Copyright (C) 2023 - 2025 Luciano Buglioni  (luciano.buglioni@gmail.com)         *
*               https://www.opensourcemech.com                                     *
*                                                                                  *
* This file is part of WeldFormFEM                                                 *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

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