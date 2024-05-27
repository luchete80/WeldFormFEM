#include "Domain_d.h"
#include "VTKWriter.h"

#include <sstream>
#include <string>
#include <iostream>

using namespace std;

namespace MetFEM{
VTKWriter::VTKWriter(Domain_d *dom, const char* fname){
  //type definition to shorten coding
	std::ostringstream oss;

	//Writing Inputs in a Log file
	string fn(fname);

	oss << "Dimension = "<< dom->m_dim << "D\n"<<endl;
  oss <<  "<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""BigEndian"">"<<endl;
  oss <<  "  <UnstructuredGrid>"<<endl;
  // FOR NOT RETURN CARRIAGE
  oss <<  "    <Piece NumberOfPoints=""" <<dom->m_node_count<< """ NumberOfCells="""<<dom->m_elem_count<<""">"<<endl;  //Note that an explicit format descriptor is needed when using
  //write (1, '(A,2x,I5)') '<Piece NumberOfPoints="'
  oss << "      <Points>"<<endl;
  oss << "        <DataArray type=""Float32"" Name=""Position"" NumberOfComponents=""3"" Format=""ascii"">"<<endl;

  
  for (int i=0;i<dom->m_elem_count;i++){
    vector_t x = dom->getPosVec(i);
    oss << x.x <<" "<<x.y <<" " <<x.z<<endl;
  }
}
};