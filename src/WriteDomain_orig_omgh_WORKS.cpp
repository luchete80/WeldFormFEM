

///////////////// ORIGINAL:WORKING
void ReMesher::WriteDomain(){
    //const int dim = m_dom->m_dim;
    
  //if (m_mesh.nelems() != m_dom->m_elem_count ||
  //    m_mesh.nverts() != m_dom->m_node_count) {
  
  
  //memcpy_t(m_->m_elnod, elnod_h, sizeof(int) * dom->m_elem_count * m_dom->m_nodxelem); 
  double *ufield  = new double [3*m_mesh.nverts()];      
  double *vfield   = new double [3*m_mesh.nverts()]; 
  double *afield   = new double [3*m_mesh.nverts()]; 
  double *pafield   = new double [3*m_mesh.nverts()];  //prev_a
  
  double *esfield  = new double [m_mesh.nelems()]; 
  double *pfield   = new double [m_mesh.nelems()]; 
  double *sigfield = new double [6*m_mesh.nelems()]; 
  double *syfield  = new double [m_mesh.nelems()]; 
  double *psfield  = new double [m_mesh.nelems()]; 
  
  
  
  double *str_rate = new double [6*m_mesh.nelems()]; 
  double *rot_rate = new double [6*m_mesh.nelems()];  
  double *tau      = new double [6*m_mesh.nelems()];   


  double *rho   = new double [m_mesh.nelems()];   

  double *vol_0 = new double [m_mesh.nelems()];    

  double *idetF   = new double [m_dom->m_elem_count];  //Inverse Deformation gradient
  
  for (int e=0;e<m_dom->m_elem_count;e++)
    idetF[e] = m_dom->vol_0[e]/m_dom->vol[e];
  
  cout << "Mapping detF"<<endl;


  
  MapElemVector<3>(m_mesh, vol_0, idetF); //Map the inverse and then multiply by the original
  cout << "vol_0 elem 0 "<<vol_0[0]<<endl;
//  cout << "MAPPING"<<endl;
  ///// BEFORE REDIMENSION!
  //MapNodalVector<3>(m_mesh, ufield,  m_dom->u);
  
  MapNodal(ufield,  m_dom->u);
  
  MapNodalVector<3>(m_mesh, vfield,  m_dom->v);
  MapNodalVector<3>(m_mesh, afield,  m_dom->a);
  MapNodalVector<3>(m_mesh, pafield,  m_dom->prev_a);
    
  MapElemVector<3>(m_mesh, esfield,   m_dom->pl_strain);
  MapElemVector<3>(m_mesh, pfield,    m_dom->p);
  //cout << "mapping press"<<endl;
  //HybridProjectionElemToElem<3>(m_mesh, pfield,    m_dom->p);

////////
/*
  //1., GET TAGS
  auto old_tag = m_mesh.get_tag<Real>(Omega_h::REGION, "plastic_strain");
  auto old_data = old_tag.array(); // Reals

  // Get the mapping from old elements to new elements
  auto elem_map = m_mesh.get_transfer_map(Omega_h::REGION); // LOs
  //3. REMAP
  Omega_h::Reals new_data = Omega_h::map_data(1, old_data, elem_map);

  //
*/  

  MapElemVector<3>(m_mesh, rho,       m_dom->rho);
  
  
  //MapElemVector<3>(m_mesh, rho_0,     m_dom->rho_0);
  double rho_0 = m_dom->rho_0[0];
  


  auto volumes = Omega_h::measure_elements_real(&m_mesh);  // One value per element
  double total_volume = Omega_h::get_sum(m_mesh.comm(), volumes);
  
  cout << "New mesh volume "<<total_volume<<endl;
  
    
  MapElemVector<3>(m_mesh, sigfield,  m_dom->m_sigma      , 6);
  MapElemVector<3>(m_mesh, str_rate,  m_dom->m_str_rate   , 6);
  MapElemVector<3>(m_mesh, rot_rate,  m_dom->m_rot_rate   , 6);
  MapElemVector<3>(m_mesh, tau,       m_dom->m_tau        , 6);

  MapElemVector<3>(m_mesh, syfield,  m_dom->sigma_y  , 1);
  MapElemVector<3>(m_mesh, psfield,  m_dom->pl_strain  , 1);
    
  /////////////////////// BOUNDARY CONDITIONS
  int bccount[3];
  int    *bcx_nod,*bcy_nod,*bcz_nod;
  double *bcx_val,*bcy_val,*bcz_val;
  for (int d=0;d<3;d++) bccount[d]=m_dom->bc_count[d];
    
  bcx_nod =new int[bccount[0]];bcx_val =new double[bccount[0]];
  bcy_nod =new int[bccount[1]];bcy_val =new double[bccount[0]]; 
  bcz_nod =new int[bccount[2]];bcz_val =new double[bccount[0]];
  
  
    
  for (int b=0;b<bccount[0];b++){bcx_nod[b]=m_dom->bcx_nod[b];bcx_val[b]=m_dom->bcx_val[b];}
  for (int b=0;b<bccount[0];b++){bcy_nod[b]=m_dom->bcy_nod[b];bcy_val[b]=m_dom->bcy_val[b];}
  for (int b=0;b<bccount[0];b++){bcz_nod[b]=m_dom->bcz_nod[b];bcz_val[b]=m_dom->bcz_val[b];}
    
  //m_dom->bc_count[0] = m_dom->bc_count[1] = m_dom->bc_count[2] = 0;
  
  
  
  ////BEFORE REWRITE
  //// WRITE
  m_dom->Free();
  
  m_dom->m_node_count = m_mesh.nverts();
  m_dom->m_elem_count = m_mesh.nelems();
  
  
  
  m_dom->SetDimension(m_dom->m_node_count,m_dom->m_elem_count);	 //AFTER CREATING DOMAIN

  malloc_t(m_dom->m_elnod, unsigned int,m_dom->m_elem_count * m_dom->m_nodxelem);

  
  cout << "SETTING DENSITY TO "<<rho_0<<endl;
  ///// TO MODIFY, SAVE AS THE MATERIAL
  m_dom->setDensity(rho_0);

  // ATENTION:
  //Volumes could also be calculated with Jacobian
  for (int e=0;e<m_dom->m_elem_count;e++){
    m_dom->vol_0[e] = vol_0[e]*volumes[e];
    //m_dom->vol  [e] = volumes[e];
    //cout << "VOL "<<e<<", "<<m_dom->vol_0[e]<<"VOLUME 0 " << volumes[e]<< ", detF" <<idetF[e]<<endl;
  }
  
  //cout << "COPYING "<<m_dom->m_elem_count * m_dom->m_nodxelem<< " element nodes "<<endl;
  memcpy_t(m_dom->u,       ufield, sizeof(double) * m_dom->m_node_count * 3);    
  memcpy_t(m_dom->v,       vfield, sizeof(double) * m_dom->m_node_count * 3);    
  memcpy_t(m_dom->a,       afield, sizeof(double) * m_dom->m_node_count * 3);   
  memcpy_t(m_dom->prev_a, pafield, sizeof(double) * m_dom->m_node_count * 3);   
    

  memcpy_t(m_dom->pl_strain, esfield,  sizeof(double) * m_dom->m_elem_count ); 
  memcpy_t(m_dom->m_sigma  ,    sigfield,   sizeof(double) * m_dom->m_elem_count *6); 
  memcpy_t(m_dom->m_str_rate,   str_rate,   sizeof(double) * m_dom->m_elem_count *6); 
  memcpy_t(m_dom->m_rot_rate,   rot_rate,   sizeof(double) * m_dom->m_elem_count *6); 
  memcpy_t(m_dom->m_tau,        tau,        sizeof(double) * m_dom->m_elem_count *6); 
  
  memcpy_t(m_dom->sigma_y,   syfield,  sizeof(double) * m_dom->m_elem_count ); 
  memcpy_t(m_dom->pl_strain, psfield,  sizeof(double) * m_dom->m_elem_count ); 

  memcpy_t(m_dom->p,          pfield,  sizeof(double) * m_dom->m_elem_count ); 
    
  memcpy_t(m_dom->rho,          rho,  sizeof(double) * m_dom->m_elem_count ); 
   
  m_dom->AssignMatAddress();
  const Material_ *matt  = &m_dom->materials[0];
  cout << "G "<<matt->Elastic().G()<<endl;
  
  
  
  double *x_h = new double [3*m_mesh.nverts()];
  int 		*elnod_h = new int [m_mesh.nelems() * 4]; //Flattened
  

  //cout << "CONVERTING MESH"<<endl;
  auto coords = m_mesh.coords();
     cout << "Setting conn"<<endl; 
  auto f = OMEGA_H_LAMBDA(LO vert) {
    auto x = get_vector<3>(coords, vert); //dim crashes
    //std::cout<< "VERT "<<vert<<std::endl;
  //for (int n = 0; n < mesh.nverts(); n++) {
    bool found = false;  // Flag to indicate whether the node is inside an element in the old mesh
    //std::cout<< "NODES "<<std::endl;
    //std::cout << m_mesh.coords()[vert]<<std::endl;
    
      // Get coordinates for the node in the new mesh
      std::array<double, 3> target_node = {x[0], x[1], x[2]}; // Now using 3D coordinates
      for (int d=0;d<3;d++)x_h[3*vert+d] = x[d];

    //n++;
  };//NODE LOOP
  parallel_for(m_mesh.nverts(), f); 

  auto elems2verts = m_mesh.ask_down(3, VERT);
        
  auto fe = OMEGA_H_LAMBDA(LO elem) {

    bool found = false;  // Flag to indicate whether the node is inside an element in the old mesh
    //std::cout<< "ELEM "<<std::endl;
    for (int ve=0;ve<4;ve++){
      auto v = elems2verts.ab2b[elem * 4 + ve];
      elnod_h[4*elem+ve] = v;
      //cout << v <<" ";
      }
    //cout <<endl;
  };//NODE LOOP  
  parallel_for(m_mesh.nelems(), fe); 
  
  cout << "NEW MESH. Done mapping "<<endl;
  cout << "Node count "<<m_dom->m_node_count<<", ELEM COUNT "<<m_dom->m_elem_count<<endl;
  memcpy_t(m_dom->x,      x_h, 3*sizeof(double) * m_dom->m_node_count);       
  memcpy_t(m_dom->m_elnod,  elnod_h, 4*sizeof(int) * m_dom->m_elem_count);  
    cout << "Setting conn"<<endl;
  m_dom->setNodElem(elnod_h); 

  
  ReMapBCs(bcx_nod,bcx_val,m_dom->bcx_nod, m_dom->bcx_val, bccount[0]);
  ReMapBCs(bcy_nod,bcy_val,m_dom->bcy_nod, m_dom->bcy_val, bccount[1]);
  ReMapBCs(bcz_nod,bcz_val,m_dom->bcz_nod, m_dom->bcz_val, bccount[2]);

  
  //cout << "deleting "<<endl;
  delete [] vfield, afield, pafield, esfield,pfield,sigfield, syfield, vol_0;
  delete [] bcx_nod,bcy_nod,bcz_nod,bcx_val,bcy_val,bcz_val;
    cout << "MESH CHANGED"<<endl;
  //} else {
      //std::cout << "Mesh is the same "<<endl;
  //}
  //cout << "Done"<<endl;
}