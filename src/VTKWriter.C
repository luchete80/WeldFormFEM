#include "Domain_d.h"
#include "VTKWriter.h"

#include <fstream>  // ofstream

using namespace std;

namespace MetFEM{
VTKWriter::VTKWriter(Domain_d *dom, const char* fname){
  //type definition to shorten coding

	//Writing Inputs in a Log file
	m_fname = fname;

	m_oss << "Dimension = "<< dom->m_dim << "D\n"<<endl;
  m_oss <<  "<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""BigEndian"">"<<endl;
  m_oss <<  "  <UnstructuredGrid>"<<endl;
  // FOR NOT RETURN CARRIAGE
  m_oss <<  "    <Piece NumberOfPoints=""" <<dom->m_node_count<< """ NumberOfCells="""<<dom->m_elem_count<<""">"<<endl;  //Note that an explicit format descriptor is needed when using
  //write (1, '(A,2x,I5)') '<Piece NumberOfPoints="'
  m_oss << "      <Points>"<<endl;
  m_oss << "        <DataArray type=""Float32"" Name=""Position"" NumberOfComponents=""3"" Format=""ascii"">"<<endl;

  
  for (int i=0;i<dom->m_node_count;i++){
    vector_t x = dom->getPosVec(i);
    m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;
    printf("         %f %f %f \n", x.x,x.y,x.z);
  }
 
  m_oss << "         </DataArray>"<<endl;
  m_oss << "       </Points>"<<endl;
  
  m_oss << "       <Cells>" <<endl;
  m_oss << "         <DataArray type=""Int32"" Name=""connectivity"" Format=""ascii"">"<<endl;

  for (int e=0;e<dom->m_elem_count;e++){
    m_oss << "         ";
    for (int en=0;en<dom->m_nodxelem;en++){
      m_oss <<dom->elnod_h[dom->m_nodxelem*e+en] <<" ";
    }
    m_oss << endl;
  }

  m_oss <<  ""; // END OF LINE
  m_oss <<  "        </DataArray>"<<endl;
  m_oss <<  "        <DataArray type=""Int32"" Name=""offsets"" Format=""ascii"">"<<endl;
  
  //TODO: CREATE A VERSION OF OFFSET
  int offs = dom->m_nodxelem;
  for (int e=0;e<dom->m_elem_count;e++){
      m_oss << offs << " ";
      //write(1,"(I10,1x)",advance="no") offs
      offs = offs + dom->m_nodxelem;
  
  }

  m_oss <<  endl; // END OF LINE
  m_oss <<  "        </DataArray>"<<endl;
  
}

void VTKWriter::writeFile(){
  string fn(m_fname);
  //fn.append("_log.dat");
	std::ofstream of(fn.c_str(), std::ios::out);
	of << m_oss.str();
	of.close();
}


};