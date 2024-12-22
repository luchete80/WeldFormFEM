#include "Domain_d.h"
#include "VTKWriter.h"

#include <fstream>  // ofstream
#include "Mesh.h"

using namespace std;

namespace MetFEM{
VTKWriter::VTKWriter(Domain_d *dom, const char* fname){
  //type definition to shorten coding

	//Writing Inputs in a Log file
	m_fname = fname;

	//m_oss << "Dimension = "<< dom->m_dim << "D\n"<<endl;
  m_oss <<  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<endl;
  m_oss <<  "  <UnstructuredGrid>"<<endl;
  // FOR NOT RETURN CARRIAGE
  int nc = dom->m_node_count;
  if (dom_d->isContactOn())
    nc += dom_d->getTriMesh()->nodecount;
  m_oss <<  "    <Piece NumberOfPoints=\"" <<nc<< "\" NumberOfCells=\""<<dom->m_elem_count<<"\">"<<endl;  //Note that an explicit format descriptor is needed when using
  //write (1, '(A,2x,I5)') '<Piece NumberOfPoints="'
  m_oss << "      <Points>"<<endl;
  m_oss << "        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" Format=\"ascii\">"<<endl;

  
  for (int i=0;i<dom->m_node_count;i++){
    if (dom->m_dim == 3){
      vector_t x = dom->getPosVec3(i);
      m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;
    } else {
      double2 x = dom->getPosVec2(i);
      m_oss << x.x <<" "<<x.y <<endl;      
      
    }
    
    //printf("         %f %f %f \n", x.x,x.y,x.z);
  }
 
  m_oss << "         </DataArray>"<<endl;
  m_oss << "       </Points>"<<endl;
  
  m_oss << "       <Cells>" <<endl;
  m_oss << "         <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">"<<endl;

  for (int e=0;e<dom->m_elem_count;e++){
    m_oss << "         ";
    for (int en=0;en<dom->m_nodxelem;en++){
      m_oss <<dom->elnod_h[dom->m_nodxelem*e+en] <<" ";
    }
    m_oss << endl;
  }

  m_oss <<  ""; // END OF LINE
  m_oss <<  "        </DataArray>"<<endl;
  m_oss <<  "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">"<<endl;
  
  //TODO: CREATE A VERSION OF OFFSET
  int offs = dom->m_nodxelem;
  for (int e=0;e<dom->m_elem_count;e++){
      m_oss << offs << " ";
      //write(1,"(I10,1x)",advance="no") offs
      offs = offs + dom->m_nodxelem;
  
  }

  m_oss <<  endl; // END OF LINE
  m_oss <<  "        </DataArray>"<<endl;


  m_oss <<  "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">"  <<endl;
  for (int e=0;e<dom->m_elem_count;e++){
    if (dom->m_dim==2){
      m_oss <<  "9 ";
    }else if (dom->m_dim==3){
      m_oss <<  "12 ";
    }
      // if (dim .eq. 2) then
        // write(1,"(I3,1x)",advance="no") 9
      // else 
        // write(1,"(I3,1x)",advance="no") 12
      // end if
  }
  m_oss <<   "" <<endl;
  m_oss <<   "        </DataArray>" <<endl;
  m_oss <<   "      </Cells>" <<endl;
  
  m_oss << "      <PointData Scalars=\"scalars\">"  <<endl;
  
  m_oss << "        <DataArray type=\"Float32\" Name=\"u\" NumberOfComponents=\""<<dom->m_dim<<"\" Format=\"ascii\"> " <<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x = dom->getDispVec(n);
    m_oss << fixed<<x.x <<" "<<x.y <<" " <<x.z<<endl;    
    
  }

  //-----
  m_oss << "    </DataArray>" <<endl;
  m_oss << "        <DataArray type=\"Float32\" Name=\"v\" NumberOfComponents=\""<<dom->m_dim<<"\" Format=\"ascii\"> " <<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x = dom->getVelVec(n);
    m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;    
    
  }
  m_oss << "    </DataArray>" <<endl;
  //------
  m_oss << "        <DataArray type=\"Float32\" Name=\"f\" NumberOfComponents=\""<<dom->m_dim<<"\" Format=\"ascii\"> " <<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x ;//= dom->getIntForceVec(n);
    m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;    
    
  }
  m_oss << "    </DataArray>" <<endl;
  //------
  
  
  m_oss << "    </PointData>" <<endl;
  
  
  m_oss << "    </Piece>" <<endl;
  m_oss << "  </UnstructuredGrid>" <<endl;
  m_oss << "</VTKFile>" <<endl;
  
  
}

void VTKWriter::writeFile(){
  string fn(m_fname);
  //fn.append("_log.dat");
	std::ofstream of(fn.c_str(), std::ios::out);
	of << m_oss.str();
	of.close();
}


};