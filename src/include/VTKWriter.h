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
