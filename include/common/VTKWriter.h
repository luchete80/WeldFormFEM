/*************************************************************************/
/*  VTKWriter.h                                                  */
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

#include <sstream>
#include <string>
#include <iostream>

////// XML STYLE
namespace MetFEM{
class VTUWriter{
public:
  VTUWriter(Domain_d *dom, const char *fname);
  void writeFile();
protected:
  std::ostringstream m_oss;
  std::string m_fname;
};


class VTKWriter{
public:
  VTKWriter(Domain_d *dom, const char *fname);
  void writeFile();
protected:
  std::ostringstream m_oss;
  std::string m_fname;
};
  
  
  
  
  
  
};
