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


#include <fstream>
#include <iomanip>	//ONY FOR GCC!!

using namespace MetFEM;

using namespace std;

#include <omp.h>


using namespace LS_Dyna;

int main(int argc, char **argv) {
 
  int m_node_count = 4;
  int m_elem_count = 1;
  int m_nodxelem   = 4;
  
  Domain_d *dom_d = new Domain_d;
  cout << "Set dimension"<<endl;
  dom_d->SetDimensionImplicit(m_node_count,m_elem_count);	 //AFTER CREATING DOMAIN
  cout << "Done. Setting nodes..."<<endl; 
  dom_d->setNode(0,0,0,0);
  dom_d->setNode(1,0.1,0,0);
  dom_d->setNode(2,0,0.1,0);
  dom_d->setNode(3,0,0,0.1);
  cout << "Done "<<endl;
  
  ////#ifdef CUDA
  // double *x_H =  new double [m_dim*m_node_count];
  // for (int n=0;n<m_node_count;n++){
    // //cout << "Node "<<n<<endl;
    // for (int d=0;d<3;d++){
      // //cout <<reader.m_node[n].m_x[d]<< " ";
      // x_H[3*n+d] = reader.m_node[n].m_x[d]; 
    // }
    // //cout <<endl;
  // }
  
  int *elnod_h       = new int [m_elem_count * m_nodxelem]; //Flattened  
  elnod_h[0]=0;elnod_h[1]=1;elnod_h[2]=2;elnod_h[3]=3;
  cout << "Setting Nod Elements"<<endl;
  dom_d->setNodElem(elnod_h);
  cout << "NODEL "<<endl;
  cout << dom_d->getElemNode(0,0)<<endl;
  
  cout << "Node Size: "<<m_node_count<<endl;  
  cout << "Element Size: "<<m_elem_count<<endl;  
  //cout << "Domain Dimension: "<< m_dim<<endl;
     
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
  for (int i=0;i<3;i++)dom_d->AddBCVelNode(0,i,0);
  dom_d->
  AddBCVelNode(1,1,0);dom_d->AddBCVelNode(1,2,0);
  dom_d->AddBCVelNode(2,0,0);dom_d->AddBCVelNode(2,2,0);
  
  //dom_d->AddBCVelNode(3,2,-1.0e-6);
  
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
  
  Solver_Eigen *solver = new Solver_Eigen();
  
  solver->setDomain(dom_d);
  solver->Allocate();
  cout << "Assemblying matrix "<<endl;
  solver->assemblyGlobalMatrix();
  cout << "Done."<<endl;
  
  solver->applyDirichletBCs();
  // cout << "Solving system"<<endl;
  solver->SetRDOF(11,-1.0e3); // EITHER FORCE OR DISP HERE
  
  //solver->Solve();

  ////// BACK TO ZERO EXTERNAL FORCE
  dom_d->SetEndTime (0.001);
  solver->SetRDOF(11,0.0);
  //dom_d->AddBCVelNode(3,2,-1.0);
  //dom_d->AllocateBCs();
  
  dom_d->SetDT(0.001); 

  
  //dom_d->setContForceVec(3,2,-1000.0);
  
  //GLOBAL MATRIX NONLINEAR METHOD
  dom_d->SolveImplicitGlobalMatrix();

	// cout << "Program ended."<<endl;
      

	
}
