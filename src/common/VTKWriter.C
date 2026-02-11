/*************************************************************************/
/*  VTKWriter.C                                                  */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                */
/*  ('https://www.opensourcemech.com',)                                    */
/*                                                                       */
/*  Copyright (c) 2023-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. */ 
/*  See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/



#include "Domain_d.h"
#include "VTKWriter.h"

#include <fstream>  // ofstream
#include "Mesh.h"
#include <iomanip>
#include "Matrix_temp.h"

#if CUDA_BUILD
#include "tensor.cuh"
#else
//#include "Matrix.h"
#endif

#include "Tensor3.C"

using namespace std;



namespace MetFEM{

void printDummyElem(Domain_d *dom, ostringstream &m_oss){ 
    if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->elemcount;n++)
        m_oss <<0.0 <<endl;       
}

void printDummyNod(Domain_d *dom, ostringstream &m_oss){
    if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
        m_oss <<0.0 <<endl;       
}

VTUWriter::VTUWriter(Domain_d *dom, const char* fname){
  //type definition to shorten coding

	//Writing Inputs in a Log file
	m_fname = fname;

	//m_oss << "Dimension = "<< dom->m_dim << "D\n"<<endl;
  cout << "Writing VTK output "<<m_fname<<endl;
  m_oss <<  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<endl;
  m_oss <<  "  <UnstructuredGrid>"<<endl;
  // FOR NOT RETURN CARRIAGE
  int nc = dom->m_node_count;
  if (dom->isContactOn())
    nc += dom->getTriMesh()->nodecount;

  int ne = dom->m_elem_count;
  if (dom->isContactOn())
    ne += dom->getTriMesh()->elemcount;
  m_oss <<  "    <Piece NumberOfPoints=\"" <<nc<< "\" NumberOfCells=\""<< ne <<"\">"<<endl;  //Note that an explicit format descriptor is needed when using
  //write (1, '(A,2x,I5)') '<Piece NumberOfPoints="'
  m_oss << "      <Points>"<<endl;
  m_oss << "        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" Format=\"ascii\">"<<endl;

  
  for (int i=0;i<dom->m_node_count;i++){
    if (dom->m_dim == 3){
      vector_t x;
      #ifdef CUDA_BUILD
      x = dom->getPosVec3_h(i);
      #else
      x = dom->getPosVec3(i);
      #endif
      m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;
    } else {
      double2 x;
      #ifdef CUDA_BUILD
      x = dom->getPosVec2_h(i);
      #else
      x = dom->getPosVec2(i);
      #endif
      m_oss << x.x <<" "<<x.y <<"0.0"<<endl;      
      
    }
    
    //printf("         %f %f %f \n", x.x,x.y,x.z);
  }
  if (dom->isContactOn()){
    for (int n=0;n<dom->getTriMesh()->nodecount;n++){
      vector_t x = dom->getTriMesh()->node[n];
      m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;      
      }
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
  
  if (dom->isContactOn()){
    for (int e=0;e<dom->getTriMesh()->elemcount;e++){
    m_oss << "         ";
    int enmax = 3;
    if (dom->m_dim == 2) enmax = 2;
      for (int en=0;en<enmax;en++){
        m_oss << dom->getTriMesh()->elnode[enmax*e+en]<<" ";
      }      
      m_oss << endl;  
    }
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
  
  if (dom->isContactOn()){
    for (int e=0;e<dom->getTriMesh()->elemcount;e++){
      m_oss << offs << " ";
      //write(1,"(I10,1x)",advance="no") offs
      offs = offs + 3;
    }
  }

  m_oss <<  endl; // END OF LINE
  m_oss <<  "        </DataArray>"<<endl;


  m_oss <<  "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">"  <<endl;
  for (int e=0;e<dom->m_elem_count;e++){
    if (dom->m_dim==2){ //TODO: CHANGE 
      if (dom->m_nodxelem == 4)
        m_oss <<  "9 ";
      else
        m_oss <<  "5 ";
    }else if (dom->m_dim==3){
      if (dom->m_nodxelem==8)
        m_oss <<  "12 ";
      else
        m_oss <<  "10 ";
    }
      // if (dim .eq. 2) then
        // write(1,"(I3,1x)",advance="no") 9
      // else 
        // write(1,"(I3,1x)",advance="no") 12
      // end if
  }
  if (dom->isContactOn()){
    for (int e=0;e<dom->getTriMesh()->elemcount;e++){
      if (dom->getTriMesh()->dimension == 3) m_oss <<  "5 ";
      else                                   m_oss <<  "3 ";
    }
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
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << fixed<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   
    
  //-----
  m_oss << "    </DataArray>" <<endl;
  m_oss << "        <DataArray type=\"Float32\" Name=\"v\" NumberOfComponents=\""<<dom->m_dim<<"\" Format=\"ascii\"> " <<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x = dom->getVelVec(n);
    m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;    
    
  }
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << fixed<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   
  m_oss << "    </DataArray>" <<endl;
  //------
  m_oss << "        <DataArray type=\"Float32\" Name=\"f\" NumberOfComponents=\""<<dom->m_dim<<"\" Format=\"ascii\"> " <<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x ;//= dom->getIntForceVec(n);
    m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;    
    
  }

  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << fixed<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   
      
  m_oss << "    </DataArray>" <<endl;
  //------
  
  
  m_oss << "    </PointData>" <<endl;
  
  
  m_oss << "    </Piece>" <<endl;
  m_oss << "  </UnstructuredGrid>" <<endl;
  m_oss << "</VTKFile>" <<endl;
  
  
}



VTKWriter::VTKWriter(Domain_d *dom, const char* fname){
  //type definition to shorten coding

	//Writing Inputs in a Log file
	m_fname = fname;

	//m_oss << "Dimension = "<< dom->m_dim << "D\n"<<endl;
  cout << "Writing VTK output "<<m_fname<<endl;
  m_oss <<  "# vtk DataFile Version 3.0"<<endl;
  m_oss <<  "VTK Example"<<endl;
  m_oss <<  "ASCII"<<endl;
  m_oss <<  "DATASET UNSTRUCTURED_GRID"<<endl;
  m_oss <<  "POINTS ";
  
  
  // FOR NOT RETURN CARRIAGE
  int nc = dom->m_node_count;
  if (dom->isContactOn()){
    nc += dom->getTriMesh()->nodecount;
  }
  m_oss << nc << " float "<<endl;
  
  int ne = dom->m_elem_count;
  if (dom->isContactOn()){
    ne += dom->getTriMesh()->elemcount;
  } 
  
  cout << "Writing nodes "<<endl;
  for (int i=0;i<dom->m_node_count;i++){
    if (dom->m_dim == 3){
      vector_t x = dom->getPosVec3_h(i);
      //for (int d=0;d<3;d++){
          m_oss << std::scientific<<std::setprecision(4)<<float(x.x) <<" "<<float(x.y) <<" " <<float(x.z)<<endl;
      //}
      //m_oss<<endll;
      //printf("Node %d %f %f %f \n", i, x.x,x.y,x.z);
    } else {
      double2 x = dom->getPosVec2_h(i);
      m_oss << x.x <<" "<<x.y <<" 0.0"<<endl;      
      
    }
    

  }
  if (dom->isContactOn()){
    for (int n=0;n<dom->getTriMesh()->nodecount;n++){
      vector_t x = dom->getTriMesh()->node[n];
      m_oss << x.x <<" "<<x.y <<" " <<x.z<<endl;      
      }
    }
 
  //TODO: MODIFY IF NOT TRIAS
  
  m_oss << "CELLS "<<ne<<" ";

  //cout << "Writing cells "<<endl;  
  int items =  (dom->m_nodxelem+1)*dom->m_elem_count ; //MODIFY WHEN NODXELEM NOT UNIFORM
  //cout << "items"<<items<<endl;
  if (dom->isContactOn()){
    //cout << "CONTACT "<<endl;
    if (dom->m_dim == 3)
      items += 4 * dom->getTriMesh()->elemcount; // ELEMENT SIDES (3 IN TRI) + 1
    else
      items += 3 * dom->getTriMesh()->elemcount; // ELEMENT SIDES (2 IN LINE) + 1
  }
  m_oss << items << endl;
  
  
  //IF REMESH 
  //delete dom->elnod_h;
  //dom->elnod_h = new int[sizeof(int) * dom->m_nodxelem * dom->m_elem_count];
  //memcpy_t(dom->elnod_h, dom->m_elnod, sizeof(int) * dom->m_nodxelem * dom->m_elem_count);   
  
  cout << "cell loop"<<endl;
  for (int e=0;e<dom->m_elem_count;e++){
    //    cout << "Element "<<e<<endl;
    m_oss << dom->m_nodxelem<<" ";
    for (int en=0;en<dom->m_nodxelem;en++){
  #ifdef CUDA_BUILD      
      m_oss <<dom->elnod_h[dom->m_nodxelem*e+en] <<" ";
  #else
      m_oss <<dom->m_elnod[dom->m_nodxelem*e+en] <<" ";    
  #endif
    }
    m_oss << endl;
  }

  int in = dom->getNodeCount();
  if (dom->isContactOn()){
    const int nen = (dom->getTriMesh()->dimension == 3) ? 3 : 2;  
    for (int e=0;e<dom->getTriMesh()->elemcount;e++){
      m_oss << nen <<" ";
      for (int en=0;en<nen;en++){
        m_oss << dom->getTriMesh()->elnode[nen*e+en]+in<<" ";
      }      
      m_oss << endl;  
    }
  }
  //cout << "Writing cell types "<<endl;
  m_oss << "CELL_TYPES "<<ne<<endl;
  for (int e=0;e<dom->m_elem_count;e++){
    //cout << "Element "<<e<<endl;
    if (dom->m_dim==2){ //TODO: CHANGE 
      if (dom->m_nodxelem == 4)
        m_oss <<  "9 ";
      else
        m_oss <<  "5 ";
    }else if (dom->m_dim==3){
      if (dom->m_nodxelem==8)
        m_oss <<  "12 ";
      else
        m_oss <<  "10 ";
    }

    m_oss <<endl;
  }
  //Contact Element Type
  if (dom->isContactOn()){
    for (int e=0;e<dom->getTriMesh()->elemcount;e++){
      if (dom->getTriMesh()->dimension == 3) m_oss <<  "5 "<<endl;  //TRI
      else                                   m_oss <<  "3 "<<endl;  //LINE
    }
  }
  

  //cout << "Writing point data "<<endl;
  double eps = 1.0e-20;
  m_oss<<"POINT_DATA "<<nc<<endl;  
  m_oss<<"VECTORS Position float"<<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x = dom->getPosVec3(n);
    if (std::abs(x.x) < eps) x.x = 0.0;
    if (std::abs(x.y) < eps) x.y = 0.0;
    if (std::abs(x.z) < eps) x.z = 0.0;
    
    if (norm(x)>1.e-10)
      m_oss << std::scientific<<std::setprecision(4)<<x.x <<" "<<x.y <<" " <<x.z<<endl;    
    else
     m_oss << std::scientific<<std::setprecision(4)<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;      
  } 
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << fixed<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   

  m_oss<<"VECTORS DISP float"<<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x = dom->getDispVec(n);
    if (norm(x)>1.e-10)
      m_oss << std::scientific<<std::setprecision(4)<<x.x <<" "<<x.y <<" " <<x.z<<endl;    
    else
     m_oss << std::scientific<<std::setprecision(4)<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;      
  } 
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << fixed<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   
  

  m_oss<<"VECTORS Acceleration float"<<endl;
  for (int n = 0; n < dom->m_node_count; n++) {
      vector_t a = dom->getAccVec(n);

      // Limpiar componentes individuales si son más chicas que eps
      if (std::abs(a.x) < eps) a.x = 0.0;
      if (std::abs(a.y) < eps) a.y = 0.0;
      if (std::abs(a.z) < eps) a.z = 0.0;
      if (dom->m_dim == 2 ) a.z = 0.0;

      // También podés aplicar length si querés ignorar todo el vector si es chico
      if (length(a) < eps)
          m_oss << "0.0 0.0 0.0" << endl;
      else
          m_oss << std::scientific << std::setprecision(4) << a.x << " " << a.y << " " << a.z << endl;
  }
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << fixed<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   

  m_oss<<"VECTORS Velocity float"<<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t v = dom->getVelVec(n);
    if (dom->m_dim==2) v.z = 0.0;
    m_oss << std::scientific<<std::setprecision(4)<<v.x <<" "<<v.y <<" " <<v.z<<endl;    
  }
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << std::scientific<<std::setprecision(4)<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   

  m_oss<<"SCALARS Part_ID float 1"<<endl;
  m_oss<<"LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_node_count;n++)
    m_oss << 0 <<endl;    
  
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << 1 <<endl;   

  if (dom->m_thermal){
    m_oss<<"SCALARS Temp float 1"<<endl;
    m_oss<<"LOOKUP_TABLE default"<<endl;
    for (int n=0;n<dom->m_node_count;n++)
      m_oss << dom->T[n] <<endl;    
    
    if (dom->isContactOn())
      for (int n=0;n<dom->getTriMesh()->nodecount;n++)
        m_oss << fixed<<dom->getTriMesh()->T_const <<endl;   
  }

  m_oss<<"VECTORS ContForce float"<<endl;
  for (int n=0;n<dom->m_node_count;n++){
    vector_t x = dom->getContForceVec(n);
    m_oss << fixed<<x.x <<" "<<x.y <<" " <<x.z<<endl;    
  }
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << fixed<<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   

  m_oss<<"SCALARS nod_mass float 1"<<endl;
  m_oss<<"LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_node_count;n++){
    m_oss << std::scientific<<std::setprecision(4)<<dom->m_mdiag[n]<<endl;    
  }
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      m_oss << 0.0<<endl;   



  
  m_oss<<"SCALARS stress float 1"<<endl;
  m_oss<<"LOOKUP_TABLE default"<<endl;
  
  
  double *a = new double[dom->m_node_count*6];
  
  //TODO: IN GPU SHOULD BE WITH HOST ARRAY INSTEAD OF DEVICE vars 
  avgScalar(dom->m_sigma,a,6)
  
  double seq;

  //Only of relevance if parts do not flow
  for (int n=0;n<dom->m_node_count;n++){

    tensor3 sig = FromFlatSym(a,n*6);
    
    tensor3 s = sig - 1.0/3.0*Trace(sig) * Identity();
    double J2 = 0.5*(s.xx*s.xx +  2.0*s.xy*s.xy + 
                                      2.0*s.xz*s.xz + 
                     s.yy*s.yy+  
                                      2.0*s.yz*s.yz +               
                     s.zz*s.zz                 
                                     
                    );
     m_oss << sqrt(3.0*J2)<<endl;
  }   
  
    printDummyNod(dom,m_oss);

  ///////////////// DEBUG
  m_oss<<"SCALARS ext_nodes float 1"<<endl;
  m_oss<<"LOOKUP_TABLE default"<<endl;
  int count = 0;
  for (int n=0;n<dom->m_node_count;n++){
    double ext = 0.0;
    if (dom->ext_nodes[n]) {ext = 1.0;count ++;}
     m_oss << ext<<endl;
  }   
  printDummyNod(dom,m_oss);

  
    m_oss<<"SCALARS nod_area float 1"<<endl;
    m_oss<<"LOOKUP_TABLE default"<<endl;
    for (int n=0;n<dom->m_node_count;n++)
      m_oss << std::scientific<<std::setprecision(6)<<dom->node_area[n] <<endl;    
    printDummyNod(dom,m_oss);

    m_oss<<"SCALARS nod_p float 1"<<endl;
    m_oss<<"LOOKUP_TABLE default"<<endl;
    for (int n=0;n<dom->m_node_count;n++)
    if (dom->p_node[n]*dom->p_node[n]>eps)
      m_oss << std::scientific<<std::setprecision(4)<<dom->p_node[n] <<endl;  
    else
      m_oss << std::scientific<<std::setprecision(4)<<"0.0"<<endl;  

        
    printDummyNod(dom,m_oss);

  m_oss<<"TENSORS SIGMAT float"<<endl;
  for (int n=0;n<dom->m_node_count;n++){
    tensor3 sig = FromFlatSym(a,n*6);
    m_oss << sig.xx << " "<<sig.xy << " "<<sig.xz<<endl;    
    m_oss << 0.0 << " "<<sig.yy << " "<<sig.yz<<endl;    
    m_oss << 0.0 << " "<<0.0 << " "<<sig.zz<<endl;        
  }
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      for (int t=0;t<3;t++)
        m_oss <<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   


  // avgScalar(dom->m_eps,a,6)
  
  // //cout << "done"<<endl;
  // m_oss<<"TENSORS EPS float"<<endl;
  // for (int n=0;n<dom->m_node_count;n++){
    // tensor3 eps = FromFlatSym(a,n*6);
    // m_oss << eps.xx << " "<<eps.xy << " "<<eps.xz<<endl;    
    // m_oss << 0.0 << " "<<eps.yy << " "<<eps.yz<<endl;    
    // m_oss << 0.0 << " "<<0.0 << " "<<eps.zz<<endl;        
  // }
  // if (dom->isContactOn())
    // for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      // for (int t=0;t<3;t++)
        // m_oss <<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   


  avgScalar(dom->m_str_rate,a,6)
  
  m_oss<<"TENSORS EPSR float"<<endl;
  for (int n=0;n<dom->m_node_count;n++){
    tensor3 eps = FromFlatSym(a,n*6);
    m_oss << eps.xx << " "<<eps.xy << " "<<eps.xz<<endl;    
    m_oss << 0.0 << " "<<eps.yy << " "<<eps.yz<<endl;    
    m_oss << 0.0 << " "<<0.0 << " "<<eps.zz<<endl;        
  }
  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->nodecount;n++)
      for (int t=0;t<3;t++)
        m_oss <<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   

  //cout <<"tensors done "<<endl;
  delete[] a;

  m_oss << "CELL_DATA "<<ne<<endl;

  m_oss<<"SCALARS ele_area float 1"<<endl;
  m_oss<<"LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss << std::scientific<<std::setprecision(6)<<dom->m_elem_area[n] <<endl;    
  printDummyElem(dom,m_oss);
  
  m_oss << "SCALARS pressure float 1"<<endl;
  m_oss << "LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss <<dom->p[n]<<endl;  

  printDummyElem(dom,m_oss);
     
  //ADD ALSO DU;;Y pressure CONTACT PARTICLES

  m_oss << "SCALARS pl_strain float 1"<<endl;
  m_oss << "LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss <<dom->pl_strain[n]<<endl;  

  printDummyElem(dom,m_oss);
  //cout << "Dummy done "<<endl;

  // m_oss << "SCALARS detJ float 1"<<endl;
  // m_oss << "LOOKUP_TABLE default"<<endl;
  // for (int n=0;n<dom->m_elem_count;n++)
    // m_oss <<dom->m_detJ[n]<<endl;  

  // printDummyElem(dom,m_oss);

  m_oss<<"TENSORS DDEVT float"<<endl;    
    for (int e=0;e<dom->m_elem_count;e++){
    tensor3 D_local= FromFlatSym(dom->m_str_rate, 6*e);

     double tr_local = Trace(D_local);
      tensor3 dev_local = D_local - Identity() * (tr_local / 3.0);
      m_oss << D_local.xx << " "<<D_local.xy << " "<<D_local.xz<<endl;    
      m_oss << 0.0 << " "<<D_local.yy << " "<<D_local.yz<<endl;    
      m_oss << 0.0 << " "<<0.0 << " "<<D_local.zz<<endl;        
    }

  if (dom->isContactOn())
    for (int n=0;n<dom->getTriMesh()->elemcount;n++)
      for (int t=0;t<3;t++)
        m_oss <<0.0 <<" "<<0.0 <<" " <<0.0<<endl;   
        
  m_oss << "SCALARS J float 1"<<endl;
  m_oss << "LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss <<dom->vol[n]/dom->vol_0[n]<<endl;  

  printDummyElem(dom,m_oss);

  m_oss << "SCALARS Vol float 1"<<endl;
  m_oss << "LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss << std::scientific<<std::setprecision(6)<<dom->vol[n]<<endl;  

  printDummyElem(dom,m_oss);

  m_oss << "SCALARS Rho float 1"<<endl;
  m_oss << "LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss <<dom->rho[n]<<endl;  

  printDummyElem(dom,m_oss);

  m_oss << "SCALARS Vol_0 float 1"<<endl;
  m_oss << "LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss<< std::scientific<<std::setprecision(6)<<dom->vol_0[n]<<endl;  

  printDummyElem(dom,m_oss);

  m_oss << "SCALARS sigy float 1"<<endl;
  m_oss << "LOOKUP_TABLE default"<<endl;
  for (int n=0;n<dom->m_elem_count;n++)
    m_oss<< std::scientific<<std::setprecision(6)<<dom->sigma_y[n]<<endl;  

  printDummyElem(dom,m_oss);
  
  cout << "DONE writing VTK."<<endl;
  
}


void VTUWriter::writeFile(){
  string fn(m_fname);
  //fn.append("_log.dat");
	std::ofstream of(fn.c_str(), std::ios::out);
	of << m_oss.str();
	of.close();
}

void VTKWriter::writeFile(){
  string fn(m_fname);
  //fn.append("_log.dat");
	std::ofstream of(fn.c_str(), std::ios::out);
	of << m_oss.str();
	of.close();
}



};
