/*************************************************************************/
/*  main.C                                                       */
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

  TriMesh_d *msh = new TriMesh_d();
  
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


  //////////////////////////////////////////////////////
  //////////////////// BOUNDARY CONDITIONS
  cout << "Setting Boundary Conditions "<<endl;
  std::vector<boundaryCondition> bConds;
  
  
    int bc_count = 0;
    //std::vector<boundaryCondition> bcondvec;
		for (auto& bc : bcs) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid,valuetype,var,ampid;
 

	////// BOUNDARY CONDITIONS
	double ampfactor = 1.0;
	
	bool free=true;
	boundaryCondition bcon;
	//bcon.type = 0;        //DEFAULT: VELOCITY
	//bcon.valueType = 0;   //DEFAULT: CONSTANT
	//bcon.value_ang = 0.0;
	readValue(bc["zoneId"], 	bcon.zoneId);
	readValue(bc["type"], 		bcon.type);
	//type 0 means velocity vc
	readValue(bc["valueType"], 	bcon.valueType);
	readVector(bc["value"], 	      bcon.value);      //Or value linear
	readVector(bc["valueAng"], 	    bcon.value_ang);  //Or Angular value

  double3 start,end;
  readValue(bc["id"], 		zoneid);
  readVector(bc["start"], 	start);
  readVector(bc["end"], 	end);

  //int bcn= dom_d->AddBCVelZone(start, end,bcon.value);
  //cout << bcn << " Nodes with Velocity " << bcon.value.x<<", "<<bcon.value.y<<", "<<bcon.value.z<<endl; 


	if (bcon.valueType == 0){//Constant

	} else {
	// if ( bcon.valueType == 1){ //Amplitude
	// readValue(bc["amplitudeId"], 		bcon.ampId);
	// readValue(bc["amplitudeFactor"], 	ampfactor); //DEFAULT NOW IS ZERO
	// bcon.ampFactor = ampfactor;
	}
	bConds.push_back(bcon);

  } //BCs
  //dom_d->AllocateBCs();
  cout << "Done "<<endl;
      
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
	
	dom_d->SearchExtNodes();
	////TODO: JOIN DIFFERENT MATRICES
  if (contact){
    cout << "Searching external nodes"<<endl;
    #ifdef CUDA_BUILD
    #else
    //dom_d->SearchExtNodes();
    #endif

    bool flipnormals = false;
    readValue(rigbodies[0]["flipNormals"],flipnormals);    
    int id = 0;
	readValue(rigbodies[0]["zoneId"],id);
    //AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
    cout <<"Creating plane mesh..."<<endl;
   //void TriMesh_d::AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens)
   cout <<"Mesh start: "<<start.x<<","<<start.y<<", "<<start.z<<endl;
   cout <<"Mesh dim: "<<dim_.x<<","<<dim_.y<<", "<<dim_.z<<endl;
   //(const int &id, const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens)
    msh->AxisPlaneMesh(0,  2, flipnormals, start , dim_,  partSide);
    dom_d->setTriMesh(msh);
	printf("Searching bcs for ZoneID %d..\n", id);
	for (int bc=0;bc<bConds.size();bc++){
		if (bConds[bc].zoneId==id){
		printf("BC Found for Zone ID: %d\n", id);
		printf("Applying Velocity %.3e %.3e %.3e\n", bConds[bc].value.x, bConds[bc].value.y, bConds[bc].value.z);
		for (int nc=0;nc<msh->nodecount;nc++)
			msh->node_v[nc]=bConds[bc].value;
		}
	}

  double penaltyfac = -1.0;
  msh->T_const = 20.;
  
  readValue(contact_[0]["fricCoeffStatic"], 	msh->mu_sta[0]); 
  readValue(contact_[0]["fricCoeffDynamic"], 	msh->mu_dyn[0]); 
  readValue(contact_[0]["heatCondCoeff"], 	  msh->heat_cond);
  readValue(contact_[0]["dieTemp"], 	  msh->T_const);
  readValue(contact_[0]["penaltyFactor"], 	penaltyfac); 
        
  printf("FRICTION COEFFS %.3e %.3e\n", msh->mu_sta[0],msh->mu_dyn[0]);
  if (penaltyfac >-1.0){
    
    dom_d->setContactPF(penaltyfac);
  }
  
	// for (int nc=0;nc<msh->nodecount;nc++)
		// msh->node_v[nc]=make_double3(0.,0.,0.);
	
    //THIS MESH AHOULD NOT BE DELETED 
	printf("-------------------\n");
  }
  
  cout <<"Done"<<endl;
  if (rigbodies.size() > 1){
    cout << "More than one Rigid Bodies found."<<endl; 
    double3 dim_,start;
    int partSide = 1;
    int id = 0;
	readValue(rigbodies[1]["zoneId"],id);
    readValue(rigbodies[1]["partSide"],partSide);
        
    readVector(rigbodies[1]["start"], 	start); 
    cout << "Start: "<<start.x<<", "<<start.y<<", "<<start.z<<endl;
    readVector(rigbodies[1]["dim"], 	dim_); 
    bool flipnormals = false;
    readValue(rigbodies[1]["flipnormals"],flipnormals);    
    
    TriMesh_d *m = new TriMesh_d();    
    m->AxisPlaneMesh(1,  2, flipnormals, start , dim_,  partSide);
    /// CONVERT TO ABSTRACT MESH
	printf("Searching bcs for ZoneID %d..\n", id);
	for (int bc=0;bc<bConds.size();bc++){
		if (bConds[bc].zoneId==id){
		printf("BC Found for Zone ID: %d\n", id);
		printf("Applying Velocity %.3e %.3e %.3e\n", bConds[bc].value.x, bConds[bc].value.y, bConds[bc].value.z);
		    m->SetMeshVel(bConds[bc].value); ////m_v
		// for (int nc=0;nc<msh->nodecount;nc++)
			// msh->node_v[nc]=bConds[bc].value;
		
		}
	}
	//m->SetMeshVel(make_double3(-10.,0.0,0.0));
    //THIS MESH AHOULD NOT BE DELETED 
	printf("-------------------\n");

    dom_d->addMeshData(*m);
    delete m;
    cout << "Done."<<endl;
    
    cout << "MESH NODE COUNT "<<msh->nodecount<<endl;
  }
  cout <<"Setting contact..."<<endl;
  //ONCE ALL MESH ARE INITIALIZED
  if (contact){
    //cout <<"Calculating plane pos"<<endl;    
      #ifdef CUDA_BUILD

    #else
      
    msh->CalcSpheres();  //NFAR Done Once if mesh is rigid
    #endif
    cout <<"Done."<<endl;
    //~ for (int e=0;e<msh->elemcount;e++)
      //~ printf("EL %d PPLANE %f\n", e, msh->pplane[e]);  
    dom_d->setContactOn();
  } //if contact
  cout << "Calulating min element size ..."<<endl;


  //dom_d->calcMinEdgeLength();
  //dx = dom_d->getMinLength();
  

  cout<< "Min length: "<<dx<<endl;
	//double dt = cflFactor * dx/(mat_cs);
  double dt = sim_time;
  //double dt = 0.800e-5;
  cout << "Time Step "<<dt<<endl;
  dom_d->SetDT(sim_time/1000.0); 
  cout << "End Time: "<<sim_time<<endl;
  dom_d->SetEndTime (sim_time);
  dom_d->setdtOut(out_time);
  
  if (remesh_interval != -1)
    dom_d->setRemeshInterval(remesh_interval); 
  //dom_d->SetEndTime (10.0*dt);
    
  
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

  // //AddBCVelNode(Node,axis,val)
  // for (int i=0;i<3;i++)dom_d->AddBCVelNode(0,i,0);
  // dom_d->
  // AddBCVelNode(1,1,0);dom_d->AddBCVelNode(1,2,0);
  // dom_d->AddBCVelNode(2,2,0);
  
  // dom_d->AddBCVelNode(3,2,-3.0e-8);
  
  // //AFTER THIS CALL
  // dom_d->AllocateBCs();
  
 
    
	cout << "Element Count "<<dom_d->getElemCount()<<endl;

    
  
  string soltype = "Dynamic";
	cout << "Reading Solver type" << endl; readValue(config["type"], 	soltype); //0: Box
 
  //~ if (soltype == "Static"){
    //~ cout << "Solver type is Static"<<endl;
    //~ dom_d->SolveStaticDisplacement();
  //~ } else if (soltype == "Dynamic") {
    //~ cout << "Solver Type is Dynamic "<<endl;
        //~ dom_d->SolveImplicitGlobalMatrix();
  //~ }
  
  //dom_d->Solve_Martins_NR();
  dom_d->Solve_Martins_Picard();
  
	cout << "Program ended."<<endl;
    
    
    

  } // if json
  
      
	} else {
    cout << "Please provide an input file."<<endl;
  }
	
}
