#include "Domain_d.h"

#include <sstream>
#include <string>
#include <iostream>

namespace MetFEM{
class VTKWriter{
public:
  VTKWriter(Domain_d *dom, const char *fname);
  void writeFile();
protected:
  std::ostringstream m_oss;
  std::string m_fname;
};
  
  
  
};