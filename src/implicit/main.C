/***********************************************************************************
* WeldFormFEM - A C++/CUDA library to simulate Mechanical Systems using            *
*               explicit Finite Element Methods                                    *
* Copyright (C) 2023 - 2025 Luciano Buglioni  (luciano.buglioni@gmail.com)         *
*               https://www.opensourcemech.com                                     *
*                                                                                  *
* This file is part of ForgeFormFEM                                                 *
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


#include <iostream>
#include "defs.h"

#include "lsdynaReader.h"
#include "nlohmann/json.hpp"
#include "Mesh.h"
#include <fstream>
#include <iomanip>	//ONY FOR GCC!!
#include "Input.h"


#include "Domain_d.h"

using namespace MetFEM;

using namespace std;

#include <omp.h>


using namespace LS_Dyna;

int main(int argc, char **argv) {
  
  int dim = 3;
	if (argc > 1){
		string inputFileName=argv[1];	
		//std::ifstream i(argv[1]);    
  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);

	Domain_d *dom_d;

  dom_d = new Domain_d;
    
  if(extension == "k") {
     
     
  lsdynaReader reader(inputFileName.c_str());


	
  double3 V = make_double3(0.0,0.0,0.0);
  double dx = 0.1;
	double3 L = make_double3(dx,dx,dx);  
	double r = 0.05;
	
	dom_d->CreateFromLSDyna(reader);
  
  } else if (extension == "json"){
    cout << "Opening JSON format file"<<endl;
    
    nlohmann::json j;

		string inputFileName=argv[1];	
		std::ifstream i(argv[1]);
		i >> j;
    
		nlohmann::json config 		= j["Configuration"];
		nlohmann::json material 	= j["Materials"];
		nlohmann::json domblock 	= j["DomainBlocks"];
		nlohmann::json domzones 	= j["DomainZones"];
		nlohmann::json amplitudes 	= j["Amplitudes"];
		nlohmann::json rigbodies 		= j["RigidBodies"];
    nlohmann::json contact_ 		= j["Contact"];
		nlohmann::json bcs 			= j["BoundaryConditions"];
		nlohmann::json ics 			= j["InitialConditions"];    

    double out_time,sim_time;
    int remesh_interval = -1;
    readValue(config["outTime"], out_time);
    readValue(config["simTime"], sim_time);
    readValue(config["reMeshStepInterval"], remesh_interval);
    

    double cflFactor = 0.3;
    readValue(config["cflFactor"], cflFactor);

    int np = 1;
    readValue(config["Nproc"], np);   
    omp_set_num_threads(np);
    //omp_set_dynamic(0);
    cout << "Number of threads used: "<< omp_get_max_threads()<<"np "<<np<<endl;
    #pragma omp parallel 
    {
      cout << "Thread "<<endl;
      cout << "Threads used: "<<omp_get_num_threads()<<", np parameter (via [Config]{Nproc})"<<endl; //INSIDE PARALLEL
    }


  /////////////-/////////////////////////////////////////////////////////////////////////////////
  // DOMAIN //
  ////////////
		string domtype = "Box";
  //Vec3_t start,L;    
  //TriMesh_d *mesh_d;
		//cout << "Reading Domain dim" << endl;  readVector(domblock[0]["dim"], 	L);
		cout << "Reading Domain type" << endl; readValue(domblock[0]["type"], 	domtype); //0: Box
    //cout << "Reading Domain mat id" << endl;  readValue(domblock[0]["matID"], 	matID); //0: Box
    //cout << "Grid Coordinate System" << endl;  readValue(domblock[0]["gridCoordSys"], 	gridCS); //0: Box
    //cout << "Slice Angle " << endl;  readValue(domblock[0]["sliceAngle"], 	slice_ang); //0: Box
    //readBoolVector(domblock[0]["sym"], 	sym); //0: Box
  cout << "Domain type: "<<domtype<<endl;
  dom_d->m_timeint_type = TimeInt::IMPLICIT;
  
  if (domtype == "File") {
        
    string filename = "";
    readValue(domblock[0]["fileName"], 	filename); 
    cout << "Reading Domain from Input file " << filename <<endl;  

    string extension = filename.substr(filename.find_last_of(".") + 1);
    cout << "Extension "<<extension<<endl;
    if(extension == "k") { 
      cout << "Creating domain from "<< filename<<endl;
      lsdynaReader reader(filename.c_str());
      dom_d->CreateFromLSDyna(reader);
    }
        //dom.ReadFromLSdyna(filename.c_str(), rho);
/*
        double totmass = 0.0;
        readValue(domblock[0]["totMass"], 	totmass); 
        
          cout << "calculating avg distance ..."<<endl;
          double avgdist = dom.getAvgMinDist();
          cout << "Avg particle distance "<<avgdist<<endl;
        
        cout <<"Setting smoothing length to "<<avgdist<<endl;
        for (int i=0;i<dom.Particles.Size();i++){
            dom.Particles[i]->h = avgdist*hfactor;
        }
        if (totmass != 0){
        for (int i=0;i<dom.Particles.Size();i++)
            dom.Particles[i]->Mass = totmass/dom.Particles.Size();
        } else 
          cout << "TOT  MASS UNKNOWN"<<endl;
*/          
  
      
    }//File



  int dim = 3;
  double tf = 5.0e-3;
	
  double3 V = make_double3(0.0,0.0,0.0);
  double dx = 0.006;
  double Lx = 0.1	;
  double Ly = 0.024;
  double Lz = 0.012;
	double3 L = make_double3(Lx,Ly,Lz);

	double r = dx/2.0;
  
  
      
      
	//dom_d->AddBoxLength(V,L,r,true);
 

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

  
  string mattype = "Bilinear";
  if      (mattype == "Bilinear")    {
    Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
    material_h  = new Material_(el);
    material_h->cs0 = sqrt(material_h->Elastic().BulkMod()/rho); //TODO: INSIDE MATERIAL 
    cout << "CS_0: "<<material_h->cs0<<endl;
    material_h->Ep = Ep;
    material_h->Material_model = BILINEAR;
    readValue(material[0]["yieldStress0"], 	material_h->sy0 );
    cout << "Material Constants, Et: "<<c[0]<<endl;
    // material_h->Material_model = BILINEAR;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Bilinear )); //
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Bilinear), cudaMemcpyHostToDevice);	

    dom_d->AssignMaterial(material_h); //Creating materials[] ptr

  } 
  dom_d->AssignMatAddress();  //set each element mat[e] address
  cout << "Done."<<endl;


      
  // //////////////////////////////////////////////////////////
  // ////////////////// RIGID BODIES //////////////////////////
  string rigbody_type;
  bool contact = false;
  if (readValue(rigbodies[0]["type"],rigbody_type))
  contact = true;
  double3 dim_,start;
  
  readVector(rigbodies[0]["start"], 	start); 
  readVector(rigbodies[0]["dim"], 	dim_); 
  bool flipnormals = false;
  readValue(rigbodies[0]["flipNormals"],flipnormals);
  int partSide = 1;
  readValue(rigbodies[0]["partSide"],partSide);

  if (contact){
    cout << "Searching external nodes"<<endl;
    #ifdef CUDA_BUILD
    #else
    dom_d->SearchExtNodes();
    #endif
    cout << "Done."<<endl;
    
    TriMesh_d *msh = new TriMesh_d();
    
    //AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
    cout <<"Creating plane mesh..."<<endl;
   //void TriMesh_d::AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens)
   cout <<"Mesh start: "<<start.x<<","<<start.y<<", "<<start.z<<endl;
   cout <<"Mesh dim: "<<dim_.x<<","<<dim_.y<<", "<<dim_.z<<endl;
    msh->AxisPlaneMesh(0,2, false, start , dim_,  partSide);

    msh->CalcSpheres();  //NFAR Done Once if mesh is rigid

    cout <<"Done"<<endl;
    msh->SetVel(make_double3(0.0,0.,-10.0));
    dom_d->setTriMesh(msh);
    
    dom_d->setContactOn();
  }
  cout << "Calulating min element size ..."<<endl;


  //dom_d->calcMinEdgeLength();
  //dx = dom_d->getMinLength();
  

  cout<< "Min length: "<<dx<<endl;
	double dt = cflFactor * dx/(mat_cs);
  //double dt = 0.800e-5;
  cout << "Time Step "<<dt<<endl;
  dom_d->SetDT(dt); 
  cout << "End Time: "<<sim_time<<endl;
  dom_d->SetEndTime (sim_time);
  dom_d->setdtOut(out_time);
  
  if (remesh_interval != -1)
    dom_d->setRemeshInterval(remesh_interval); 
  //dom_d->SetEndTime (10.0*dt);
  
  //////////////////// BOUNDARY CONDITIONS
    int bc_count = 0;
    //std::vector<boundaryCondition> bcondvec;
		for (auto& bc : bcs) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid,valuetype,var,ampid;
   
		}//Boundary Conditions
  
  
  
  int fixcount =0;
  int velcount =0;
  for (int i=0;i<dom_d->getNodeCount();i++){

    // if (dom_d->getPosVec3(i).z <0.0002) {
      // for (int d=0;d<3;d++)dom_d->AddBCVelNode(i,d,0);
      // fixcount++;
      // //cout << "node "<< i<<" fixed "<<endl;
    // }

    

    // //#ifdef CUDA_BUILD
    // //#else    
    // if (dom_d->getPosVec3_h(i).z > 0.03-0.0003 ) {
    // //if (dom_d->getNodePos3(i).z > 0.616-0.025 ) {
      // //dom_d->AddBCVelNode(i,0,-0.0);
      // //dom_d->AddBCVelNode(i,1,-0.0);
      // //dom_d->AddBCVelNode(i,2,-1.0e-3);
      // //cout << "Node "<<i <<" vel "<<endl;
      // velcount++;
    // }     
    // //#endif
    
  }
  //initElemArrayCPU (this,sigma_y,1,300.0e6)  

  //AddBCVelNode(Node,axis,val)
  for (int i=0;i<3;i++)dom_d->AddBCVelNode(0,i,0);
  dom_d->
  AddBCVelNode(1,1,0);dom_d->AddBCVelNode(1,2,0);
  dom_d->AddBCVelNode(2,2,0);
  // Factorizing
  // Solution:
       // 0
       // 0
       // 0
   // 9e-09
       // 0
       // 0
       // 0
   // 9e-09
       // 0
       // 0
       // 0
  // -3e-08
  dom_d->AddBCVelNode(3,2,-3.0e-8);
  
  //AFTER THIS CALL
  dom_d->AllocateBCs();
  
 
  //
  cout << "Reading contact "<<endl;

    
	cout << "Element Count "<<dom_d->getElemCount()<<endl;

    
    
    
    

  }
  
  ////// ELASTIC TEST
  //dom_d->Solve();
  dom_d->ElasticIncSolve();
	//dom_d->Solve ();
  
  
	cout << "Program ended."<<endl;
      
	} else {
    cout << "Please provide an input file."<<endl;
  }
	
}
