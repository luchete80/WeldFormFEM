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
#include "Solver_Eigen.h"
#include <iostream>
#include "defs.h"
#include "VTKWriter.h"

#include <fstream>
#include <iomanip>	//ONY FOR GCC!!

using namespace MetFEM;

using namespace std;

#include <omp.h>


using namespace LS_Dyna;
using namespace MetFEM;

int Domain::SolveElastic() {
 
  int m_node_count = 4;
  int m_elem_count = 1;
  int m_nodxelem   = 4;
  
  Domain_d *dom_d = new Domain_d;
  dom_d->m_timeint_type = TimeInt::IMPLICIT;
  
  //~ cout << "Set dimension"<<endl;
  //~ dom_d->SetDimensionImplicit(m_node_count,m_elem_count);	 //AFTER CREATING DOMAIN
  //~ cout << "Done. Setting nodes..."<<endl; 
  //~ dom_d->setNode(0,0,0,0);
  //~ dom_d->setNode(1,0.1,0,0);
  //~ dom_d->setNode(2,0,0.1,0);
  //~ dom_d->setNode(3,0,0,0.1);
  
  std::string filename = "tetra_cyl.k";
  
  lsdynaReader reader(filename.c_str());
  dom_d->CreateFromLSDyna(reader);
      
      
  cout << "Done "<<endl;
  
     
  ////// MATERIAL  
  double E, nu, rho;
  E   = 200.0e9;
  nu  = 0.3;
  rho = 7850.0;

  cout << "Setting density"<<endl;
  dom_d->setDensity(rho); //rho_0
  cout <<"Done."<<endl;
  cout << "Creating Material..:"<<endl;
  Material_ *mat_h = (Material_ *)malloc(dom_d->getElemCount() * sizeof(Material_ *)); 
  Elastic_ el(E,nu);
  // cout << "Mat type  "<<mattype<<endl;

  Material_ *material_h;
  double Ep, c[6];

  
  double mat_modK= E / ( 3.0*(1.0 -2.0*nu) );
  double mat_G= E / (2.0* (1.0 + nu));
  
  // dom%mat_K = mat_modK !!!TODO CREATE MATERIAL
  
  double mat_cs = sqrt(mat_modK/rho);

  Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
  material_h  = new Material_(el);
  material_h->cs0 = sqrt(material_h->Elastic().BulkMod()/rho); //TODO: INSIDE MATERIAL 
  cout << "CS_0: "<<material_h->cs0<<endl;
  material_h->Ep = Ep;
  material_h->Material_model = BILINEAR;
  

  dom_d->AssignMaterial(material_h);

  cout << "Done."<<endl;
  
  int fixcount =0;
  int velcount =0;
    
  //AddBCVelNode(Node,axis,val)
  for (int i=0;i<dom_d->getNodeCount();i++){
    // if (dom_d->getPosVec3_h(i).z<0.029){
      // for (int d=0;d<2;d++) dom_d->AddBCVelNode(i,d,0.0);
      // for (int d=0;d<2;d++) dom_d->AddBCVelNode(i,d,-1.0e-4);
      // velcount++;
    // }
    
    if (dom_d->getPosVec3_h(i).z<0.0005){
      for (int d=0;d<3;d++) dom_d->AddBCVelNode(i,d,0.0);
      fixcount++;
    }
  }
  //dom_d->AddBCVelNode(3,2,-0.001);
  
      //cout << "node "<< i<<" fixed "<<endl;
    
  //if (dom_d->getPosVec3_h(i).z > 0.616-0.002 ) 

  cout << "FIXED "<<fixcount<< " NODES"<<endl;  
  cout << "VEL  "<<velcount<< " NODES"<<endl;  
  
  //AFTER THIS CALL
  dom_d->AllocateBCs();

  //THIS IS LIKE SOLVE
  dom_d->AssignMatAddress();

  dom_d->calcElemJAndDerivatives();
  dom_d->CalcElemVol();

  //IMPLICIT 
  dom_d->CalcMaterialStiffElementMatrix();
  
  ////////// ALTERNATIVE, TO ASSEMBLE ALSO LARGE MATRICES
    // /////////////////////// THIS IS BEB
    // par_loop(e,m_elem_count){
          // //cout << "Element "<<e<<endl;
          
          // // 6) Build B matrix (strain-displacement) for the element
          // int tid = omp_get_thread_num();

          // Matrix &B = Bmat_per_thread[tid];
          // //// HERE B is in fact BxdetJ
          // B = getElemBMatrix(e); // dimensions 6 x (m_nodxelem * m_dim)
          // B = B *(1.0/m_detJ[e]);

          // //cout << "B mat "<<endl;
          // //B.Print();

          // // 7) Compute internal force: fint = V_e * B^T * σ

          // Matrix stress_voigt = FlatSymToVoigt(m_sigma,m_dim,m_nodxelem);
          // //CHANGE TO FORCES TO MATRIX! already calculated
          // Matrix fint = MatMul(B.getTranspose(), stress_voigt); //DO NOT TRANSPOSE B DEFITELY
 
          // fint = fint * vol[e];

          // ///TANGENT!
          // // // 8.1) Compute tangent stiffness matrix Ktan = V_e * B^T * D * B
          // Matrix D(6,6);
          // D =  mat[e]->getElasticMatrix();
          
          // Matrix Kmat = MatMul(B.getTranspose(), MatMul(D, B));
          // Kmat = Kmat * (1.0/6.0*m_detJ[e]); // B is B x detJ

          // double Ve = vol[e]; // Current volume (updated Lagrangian)

          // ///////////////////////////////////////////////////
          // /////////// IMPORTANT!!! --A LOT-- FASTER (LESS PRODUCTS) THAN: Kgeo = G^T sigma G
          // // 2. Initialize Kgeo (12x12 for 4-node tetrahedron)
          // //Matrix& Kgeo = *(m_Kgeo[e]);
          // Matrix Kgeo(m_dim*m_nodxelem,m_dim*m_nodxelem);
          // Kgeo.SetZero();
          
          // // // 3. Loop over node pairs (a, b)
          // // // REMEMBER DERIVATIVES ARE AFFECTED BY DETJ
          // for (int a = 0; a < 4; ++a) {
            // // ∇Nᵃ in current config (∂Nᵃ/∂x, ∂Nᵃ/∂y, ∂Nᵃ/∂z)
            // Matrix grad_a(3, 1);
            // grad_a.Set(0, 0, getDerivative(e, 0, 0, a)); // ∂N/∂x
            // grad_a.Set(1, 0, getDerivative(e, 0, 1, a)); // ∂N/∂y
            // grad_a.Set(2, 0, getDerivative(e, 0, 2, a)); // ∂N/∂z

            // for (int b = 0; b < 4; ++b) {
              // // ∇Nᵇ in current config
              // Matrix grad_b(3, 1);
              // grad_b.Set(0, 0, getDerivative(e, 0, 0, b));
              // grad_b.Set(1, 0, getDerivative(e, 0, 1, b));
              // grad_b.Set(2, 0, getDerivative(e, 0, 2, b));

              // // Compute K_geo(a,b) = (∇Nᵃ)ᵀ · σ · ∇Nᵇ * Ve
              // Matrix sigma_grad_b = MatMul(FlatSymToMatrix(m_sigma), grad_b); // σ · ∇Nᵇ (3x1)
              // Matrix kab = MatMul(grad_a.getTranspose(), sigma_grad_b); // 1x1 scalar
              // double k_ab = kab.getVal(0, 0) * Ve;

              // // Fill 3x3 block (assumes 3 DOF per node)
              // for (int i = 0; i < 3; ++i) {
                // Kgeo.Set(3*a + i, 3*b + i, Kgeo.getVal(3*a + i, 3*b + i) + k_ab);
              // }
            // }
          // }

          // Kgeo = Kgeo * (1.0/(6.0*m_detJ[e]));
          // Matrix K = Kgeo + Kmat;

          // K = K*dt;
          
          // double beta = 0.25;
          // // // Add mass scaling for stability (FORGE does this)
          // for (int i = 0; i < m_nodxelem; i++) {  // Loop over element nodes
              // int node = getElemNode(e, i);        // Get global node number
              // for (int d = 0; d < m_dim; d++) {   // Loop over dimensions (x,y,z)
                  // int idx = i*m_dim + d;           // Local DOF index
                  // double mass_term = m_mdiag[node] / (beta * dt);  //kg/s = (kgxm/s2) x s/m = N/m x s
                  // K.Set(idx, idx, (K.getVal(idx, idx) + mass_term) *(1.0 + 1.0e-8) ); //ALSO ADDED DIAG REGULARIZATION
              // }
          // }
          // //cout <<"CHECKING INTERNAL FORCES"<<endl;

          // Matrix R(m_dim*m_nodxelem,1);
          // for (int i = 0; i < m_nodxelem; i++) {
            // //int node = getElemNode(e, i % m_nodxelem);
            // for (int d=0;d<m_dim;d++){
            // //cout << "NODE, DIM "<<i<<","<<d<<", fint mat"<<fint.getVal(m_dim*i+d,0)<<", fel "<<m_f_elem[i*m_dim+d]<<endl;
            // //R.Set(i,0,-fint.getVal(m_dim*i+d,0)); //ADD EXTERNAL ELEMENT FORCES
            // R.Set(m_dim*i+d,0,-m_f_elem[i*m_dim+d]/*+m_f_elem_hg [offset + i*m_dim + d]*/); //ADD EXTERNAL ELEMENT FORCES
            // }
          // }
          // solver->assembleElement(e, K);
          // solver->assembleResidual(e,R);//SHOULD BE NEGATIVE!  
          
      // }//elem


  
  Solver_Eigen *solver = new Solver_Eigen();
  
  solver->setDomain(dom_d);
  solver->Allocate();
  cout << "Assemblying matrix "<<endl;
  solver->assemblyGlobalMatrix();
  cout << "Done."<<endl;
  
  solver->applyDirichletBCs();
  for (int i=0;i<dom_d->getNodeCount();i++){
      if (dom_d->getPosVec3_h(i).z>0.029){
        solver->SetRDOF(3*i+2,-10.0);
      velcount++;
    }
  }
  cout << "Applied "<<velcount <<" forces "<<endl;
  // cout << "Solving system"<<endl;
  //solver->SetRDOF(11,-1.0e3); // NOT WORKING
  
  solver->Solve();

  for (int i=0;i<dom_d->getNodeCount();i++){
    for (int d=0;d<3;d++)
      cout <<solver->getU(i,d)<<", ";
    cout <<endl;
  }

      

  VTKWriter writer(dom_d, "out.vtk");
  writer.writeFile();
	
}
