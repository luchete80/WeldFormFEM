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

#include <iostream>
#include "defs.h"

#include "lsdynaReader.h"
#include "nlohmann/json.hpp"
#include "Mesh.h"
#include <fstream>
#include <iomanip>	//ONY FOR GCC!!
#include "Input.h"

#include "git_commit.h"

using namespace MetFEM;

using namespace std;
#ifdef CUDA_BUILD
void report_gpu_mem()
{
  size_t free, total;
  cudaMemGetInfo(&free, &total);
  std::cout << "Free = " << free << " Total = " << total <<std::endl;
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#else
#include <omp.h>

#endif

size_t findLastOccurrence(string str, char ch)
{
 
    // To store the index of the result
    size_t found;

    found = str.rfind(ch);
    // If string doesn't have
    // character ch present in it
    if (found == string::npos) {
        cout << "Character " << ch
             << " is not present in"
             << " the given string.";
    }
 
    // Else print the position
    else {
        cout << "The last occurrence of '"
             << ch << "' is found at index: "
             << found << endl;
    }
  return found;
}

// Función para cargar los parámetros desde JSON
StabilizationParams loadStabilizationParams(const nlohmann::json& j,Domain_d *dom) {
    StabilizationParams params;
    
    if (j.contains("Stabilization") && !j["Stabilization"].is_null()) {
      
    // Accede al objeto "Stabilization" en el JSON
    const nlohmann::json& stab = j["Stabilization"];

    params.alpha_free       = stab.value("alpha_free", 0.0);         
    params.alpha_contact    = stab.value("alpha_contact", 0.0);      
    params.hg_coeff_free    = stab.value("hg_coeff_free", 0.0);     
    params.hg_coeff_contact = stab.value("hg_coeff_contact", 0.0);  
    params.av_coeff_div     = stab.value("av_coeff_div",  0.0);     
    params.av_coeff_bulk    = stab.value("av_coeff_bulk", 0.0);    
    params.log_factor       = stab.value("log_factor", 0.0);        
    params.p_pspg_bulkfac   = stab.value("p_pspg_bulkfac", 0.0);     
    params.J_min            = stab.value("J_min", 0.0);    
    params.hg_forces        = stab.value("hg_forces", 0.0);            
        //~ params.alpha_free       = stab.value("alpha_free", 0.0);         
    //~ params.alpha_contact    = stab.value("alpha_contact", 0.0);      
    //~ params.hg_coeff_free    = stab.value("hg_coeff_free", 0.05);     
    //~ params.hg_coeff_contact = stab.value("hg_coeff_contact", 0.05);  
    //~ params.av_coeff_div     = stab.value("av_coeff_div",  0.15);     
    //~ params.av_coeff_bulk    = stab.value("av_coeff_bulk", 0.15);    
    //~ params.log_factor       = stab.value("log_factor", 0.8);        
    //~ params.p_pspg_bulkfac   = stab.value("p_pspg_bulkfac", 0.05);    
    
    dom->m_stab = params;
    } else {
        std::cerr << "[WARNING] 'Stabilization' block not found in JSON, using defaults." << std::endl;
        // Defaults already set by constructor or initializer list
    }

    return params;
        
}
 

using namespace LS_Dyna;

int main(int argc, char **argv) {

  TriMesh_d *msh = new TriMesh_d();
  
  int dim = 3;
	if (argc > 1){
		string inputFileName=argv[1];	
		//std::ifstream i(argv[1]);    
  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);

	Domain_d *dom_d;
  dom_d->m_timeint_type == TimeInt::EXPLICIT;
  #ifdef CUDA_BUILD
	report_gpu_mem();
	gpuErrchk(cudaMallocManaged(&dom_d, sizeof(MetFEM::Domain_d)) );
	report_gpu_mem();
  #else
    dom_d = new Domain_d;
  #endif    
  

  if(extension == "k") {
     
     
  lsdynaReader reader(inputFileName.c_str());

  //TRIMESH DATA. IF MESH IS NOT RESIZED, THIS IS THE ADRESS (DO NOT DELETE)  

	
  double3 V = make_double3(0.0,0.0,0.0);
  double dx = 0.1;
	double3 L = make_double3(dx,dx,dx);  
	double r = 0.05;
	
	dom_d->CreateFromLSDyna(reader);
  
  } else if (extension == "json"){

    size_t pos = 0;
    size_t test = findLastOccurrence(inputFileName, '\\');
    if (test != string::npos) pos = test;
    
    string out_name = inputFileName.substr(pos, inputFileName.find(".json") - pos +1) + "out";
    cout << "Out file: "<< out_name << endl;
    //dom_d->out_file.open(out_name.c_str(), std::ofstream::out | std::ofstream::app);
    dom_d->out_file.open(out_name.c_str());

    ostringstream oss_out;
    
    oss_out << "-----------------------------------------------"<<endl;
    oss_out << "------------------ WELDFORM FEM -----------------"<<endl;
    oss_out << "------------- Version: " <<  PROJECT_VERSION << "------------------"<<endl;
    oss_out << "------------- Date:  " << GIT_COMMIT_HASH << "-------------------"<<endl;
    oss_out << "----------------- "<< BUILD_DATE << "-------------------"<<endl;
            
    cout << oss_out.str();
    dom_d->out_file << oss_out.str();    
    dom_d->out_file.flush(); // Optional

    cout << "Opening JSON format file"<<endl;
    
    nlohmann::json j;

		string inputFileName=argv[1];	
		std::ifstream i(argv[1]);
		i >> j;
    
		nlohmann::json config 		= j["Configuration"];
		nlohmann::json stab 		=  j["Stabilization"];
		nlohmann::json material 	= j["Materials"];
		nlohmann::json domblock 	= j["DomainBlocks"];
		nlohmann::json domzones 	= j["DomainZones"];
		nlohmann::json amplitudes 	= j["Amplitudes"];
		nlohmann::json rigbodies 		= j["RigidBodies"];
    nlohmann::json contact_ 		= j["Contact"];
		nlohmann::json bcs 			= j["BoundaryConditions"];
		nlohmann::json ics 			= j["InitialConditions"];   
		nlohmann::json mesh 			= j["Meshing"];   
    
    loadStabilizationParams(j, dom_d); 

    double out_time,sim_time;
    bool fixedTS = false;
    
    int remesh_interval = -1;
    readValue(config["outTime"], out_time);
    readValue(config["simTime"], sim_time);

    readValue(config["plHeatFrac"], dom_d->m_plheatfraction);
    
    dom_d->out_file << "Plastic Heat Fraction Coeff: " <<dom_d->m_plheatfraction<<endl;
    dom_d->out_file.flush(); // Optional
    
        
    readValue(config["fixedTS"], fixedTS);
    
    double cflFactor = 0.3;
    readValue(config["cflFactor"], cflFactor);
    dom_d->setCFL(cflFactor);
    if (fixedTS) dom_d->setFixedDt(true);


    readValue(mesh["stepInterval"], remesh_interval);
    cout << "MESH INTERVAL: "<<remesh_interval<<endl;

    readValue(mesh["minStrain"], dom_d->m_remesh_min_pl_strain);
    readValue(mesh["maxStrain"], dom_d->m_remesh_max_pl_strain);
    readValue(mesh["mapVel"],     dom_d->m_remesh_map_vel);
    readValue(mesh["mapAcc"],     dom_d->m_remesh_map_acc);
    readValue(mesh["maxCount"],   dom_d->m_remesh_max_count);
    readValue(mesh["dampFactor"],   dom_d->m_remesh_damp_vel);
    
    double3 artifvisc;
    dom_d->m_artifvisc[0] = 0.0; //Alpha. linear 0.3-0.6
    dom_d->m_artifvisc[1] = 0.0;//Beta, quadratic 0.02-0.05
    
    readVector(config["artifViscCoeffs"], artifvisc);      //Or value linear
    dom_d->m_artifvisc[0] = artifvisc.x;dom_d->m_artifvisc[1] = artifvisc.y;
    
    string dom_type = "3D";
    readValue(config["domType"], 	dom_type); 
    bool xyzsym[] = {false,false,false};
    double symtol = 1.0e-4;
    readValue(config["symtol"],symtol);
    readValue(config["xSymm"], 	xyzsym[0]);
    readValue(config["ySymm"], 	xyzsym[1]); 
    readValue(config["zSymm"], 	xyzsym[2]); 
    string symaxis[]={"X","Y","Z"};
    for (int d=0;d<3;d++)
      if (xyzsym[d])
        cout << "SYMMETRY ON AXIS "<<symaxis[d]<<endl;
     
    if (dom_type == "AxiSymm" || dom_type == "AxiSym"){
      dom_d->setAxiSymm();
      cout << "DOMAIN TYPE: AXIS SYMMETRIC"<<endl;
    } else if (dom_type == "plStrain"){
      cout << "DOMAIN TYPE: PLAIN STRAIN"<<endl;
    }
    
    int press_alg = 0;
    readValue(config["pressAlgorithm"], 	press_alg); 
    if (press_alg>0)
      dom_d->m_press_algorithm = press_alg;
    
    dom_d->out_file<<"PRESSURE ALGORITHM IS ";
    if (press_alg == 0)dom_d->out_file<<"Hybrid-FBAR"<<endl;
    if (press_alg == 1)dom_d->out_file<<"Bonet  FBAR"<<endl;
    
    dom_d->out_file.flush();
    
    
    #ifdef CUDA_BUILD
    
    #else
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
    #endif

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
    else if (domtype == "Box"){
       cout << "Adding Box ..."<<endl;  
       double3 start,L;
       double dx = 0.06;
			readVector(domblock[0]["dim"], 		L);
			readVector(domblock[0]["start"], 	start);
      readValue(domblock[0]["elemLength"], 	dx);
      bool tritet = false;
      string eltype = "";
      readValue(domblock[0]["elemType"], 	eltype);
      if (eltype == "TriTet"){
        tritet = true;
        cout << "Element type set to TRI/TET"<<endl;
      }
      cout << "Box Start: "<<start.x<< ", "<<start.y<< ", "<<start.z<<endl;
      cout << "Box Length : "<<start.x<< ", "<<start.y<< ", "<<start.z<<endl;

    // Domain_d::AddBoxLength(vector_t const & V, vector_t const & L, const double &r,const bool &red_int, const bool &tritetra)
    #ifdef CUDA_BUILD
      cout << "STILL NOT AVAIABLE"<<endl;
    #else
      dom_d->AddBoxLength(start, make_double3(L.x,L.y,0.0), dx/2.,true,tritet);	
    #endif
    }

		double IniTemp = 20.;
    int ics_count = 0;
		for (auto& ic : ics){
			double temp;
      int id;
      readValue(ic["Temp"], IniTemp);
      readValue(ic["id"], id);
      cout << "Initial Temp: "<<IniTemp<<endl;      
    }
        
  bool thermal = false;
  readValue(config["thermal"], thermal);
  if (thermal){
    dom_d->setThermalOn();
    dom_d->setTemp(IniTemp);
  }
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
  cout << "Density.."<< endl; readValue(material[0]["density0"], 		rho);
  readValue(material[0]["youngsModulus"], 	E);
  readValue(material[0]["poissonsRatio"], 	nu);

  std::vector<double> e_range (2,0.0);
  std::vector<double> er_range(2,0.0);
  std::vector<double> T_range (2,0.0);
  e_range[1]=er_range[1]=T_range[1]=1.0e10;  

  cout << "Setting density"<<endl;
  dom_d->setDensity(rho); //rho_0
  cout <<"Done."<<endl;
  cout << "Creating Material..:"<<endl;
  Material_ *mat_h = (Material_ *)malloc(dom_d->getElemCount() * sizeof(Material_ *)); 
  Elastic_ el(E,nu);
  // cout << "Mat type  "<<mattype<<endl;

  Material_ *material_h = nullptr;
  double Ep;
  std::vector<double> c;
  c.resize(10);
  
  double mat_modK= E / ( 3.0*(1.0 -2.0*nu) );
  double mat_G= E / (2.0* (1.0 + nu));
  
  // dom%mat_K = mat_modK !!!TODO CREATE MATERIAL
  
  double mat_cs = sqrt(mat_modK/rho);

  readArray(material[0]["strRange"],  e_range );
  readArray(material[0]["strdotRange"], er_range);
  readArray(material[0]["tempRange"],   T_range );  

      
  string mattype = "Bilinear";
  cout << "Type: "; 
  readValue(material[0]["type"], 		mattype);
  cout << mattype<<endl;
  double Fy=0.0;
  readValue(material[0]["yieldStress0"],Fy);
  readArray(material[0]["const"], 		c);
      
  if      (mattype == "Bilinear")    {
    Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
    material_h  = new Material_(el);
    material_h->cs0 = sqrt(material_h->Elastic().BulkMod()/rho); //TODO: INSIDE MATERIAL 
    cout << "CS_0: "<<material_h->cs0<<endl;
    material_h->Ep = Ep;
    material_h->Material_model = BILINEAR;

    // cout << "Material Constants, Et: "<<c[0]<<endl;
    // material_h->Material_model = BILINEAR;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Bilinear )); //
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Bilinear), cudaMemcpyHostToDevice);	

  } 


   else if (mattype == "Hollomon")    {
    material_h  = new Hollomon(el,Fy,c[0],c[1]);
    cout << "Hollomon Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
    // // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Hollomon));
    
    // material_h  = new Material_(el);
    material_h->InitHollomon(el,Fy,c[0],c[1]);
    material_h->Material_model = HOLLOMON;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Material_));
    
    // //init_hollomon_mat_kernel<<<1,1>>>(dom_d); //CRASH
    // //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Hollomon*), cudaMemcpyHostToDevice);	
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Material_), cudaMemcpyHostToDevice);	 //OR sizeof(Hollomon)??? i.e. derived class
    
  
  } else if (mattype == "JohnsonCook") {
    // //Order is 
                               // //A(sy0) ,B,  ,C,   m   ,n   ,eps_0,T_m, T_transition
   material_h  = new JohnsonCook(el,Fy, c[0],c[1],c[3],c[2],c[6], c[4],c[5]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
   cout << "Johnson Cook Material"<<endl; 
   
    // //Only 1 material to begin with
    // //cudaMalloc((void**)&dom_d->materials, 1 * sizeof(JohnsonCook ));
    // //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(JohnsonCook), cudaMemcpyHostToDevice);	
    // cout << "Material Constants, B: "<<c[0]<<", C: "<<c[1]<<", n: "<<c[2]<<", m: "<<c[3]<<", T_m: "<<c[4]<<", T_t: "<<c[5]<<", eps_0: "<<c[6]<<endl;
  } else if (mattype == "GMT") {
      //Order of input is: n1,n2  c1,c2, m1,m2, I1, I2
      material_h = new GMT(el,c[0],c[1],c[2],c[3],c[4], c[5],c[6],c[7],
                       e_range [0],e_range [1],
                       er_range[0],er_range[1],
                       T_range [0],T_range [1]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
      cout << "GMT Material Constants: "<<endl<<
                                  "n1: "<<c[0]<<", n2: "<<c[1]<<endl<<
                                  "c1: "<<c[2]<<", c2: "<<c[3]<<endl<<
                                  "m1: "<<c[4]<<", m2: "<<c[5]<<endl<<
                                  "I1: "<<c[6]<<", I2: "<<c[7]<<endl;

      cout << "Strain Range: "      <<material_h->e_min<<", "<<material_h->e_max <<endl;
      cout << "Strain Rate Range: " <<material_h->er_min<<", "<<material_h->er_max <<endl;
      cout << "Temp   Range: "      <<material_h->T_min<<", "<<material_h->T_max <<endl;

    }    
    
    else                              //throw new Fatal("Invalid material type.");                            
    printf("ERROR: Invalid material type.\n");

  if (material_h){
    ////// THERMAL 
    material_h->cs0 = sqrt(material_h->Elastic().BulkMod()/rho); //TODO: INSIDE MATERIAL 
    readValue(material[0]["thermalCond"], 	  material_h->k_T);
    readValue(material[0]["thermalHeatCap"], 	material_h->cp_T);    
    readValue(Fy, 	material_h->sy0 );
    dom_d->AssignMaterial(material_h);
  }
  cout << "Done."<<endl;
  
  //////////////////////////////////////////////////////
  //////////////////// BOUNDARY CONDITIONS
  std::vector<boundaryCondition> bConds;
  
  cout << "Reading Boundary Conditions ..."<<endl;
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
      
      int bcn= dom_d->AddBCVelZone(start, end,bcon.value);
      cout << bcn << " Nodes with Velocity " << bcon.value.x<<", "<<bcon.value.y<<", "<<bcon.value.z<<endl; 
  
	if (bcon.valueType == 0){//Constant

	} else {
	// if ( bcon.valueType == 1){ //Amplitude
	// readValue(bc["amplitudeId"], 		bcon.ampId);
	// readValue(bc["amplitudeFactor"], 	ampfactor); //DEFAULT NOW IS ZERO
	// bcon.ampFactor = ampfactor;
	}
	bConds.push_back(bcon);

  } //BCs
  
      
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
  readValue(rigbodies[0]["flipnormals"],flipnormals);
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
    readValue(rigbodies[0]["flipnormals"],flipnormals);    
    int id = 0;
	readValue(rigbodies[0]["zoneId"],id);
    //AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
    cout <<"Creating plane mesh..."<<endl;
   //void TriMesh_d::AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens)
   cout <<"Mesh start: "<<start.x<<","<<start.y<<", "<<start.z<<endl;
   cout <<"Mesh dim: "<<dim_.x<<","<<dim_.y<<", "<<dim_.z<<endl;
   //(const int &id, const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens)
    string rigbody_type;
    bool contact = false;
    readValue(rigbodies[0]["type"],rigbody_type);

      
    if (rigbody_type == "Plane"){
      msh->AxisPlaneMesh(0,  2, flipnormals, start , dim_,  partSide);
      dom_d->setTriMesh(msh);
    }
    else if (rigbody_type == "File"){
      string filename = "";
      readValue(rigbodies[0]["fileName"], 	filename); 
      NastranReader reader(filename.c_str());
      *msh = TriMesh_d(reader,flipnormals);
      //mesh.push_back (new SPH::TriMesh(reader,flipnormals ));  
      
    }
    
  printf("Searching bcs for ZoneID %d..\n", id);
	for (int bc=0;bc<bConds.size();bc++){
		if (bConds[bc].zoneId==id){
		printf("BC Found for Zone ID: %d\n", id);
		printf("Applying Velocity %.3e %.3e %.3e\n", bConds[bc].value.x, bConds[bc].value.y, bConds[bc].value.z);
		for (int nc=0;nc<msh->nodecount;nc++)
			msh->node_v[nc]=bConds[bc].value;
		}
	}
  cout << "NODE COUNT "<<msh->nodecount<<endl;
  
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
	
  //// TODO: CHANGE THIS TO GENERAL CONTACT
    //THIS MESH SHOULD NOT BE DELETED 
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
    
    string rigbody_type;
    bool contact = false;
    readValue(rigbodies[1]["type"],rigbody_type);
    
    if (rigbody_type == "Plane") { 
      
      m->AxisPlaneMesh(1,  2, flipnormals, start , dim_,  partSide);

    }
    else if (rigbody_type == "File"){
      string filename = "";
      readValue(rigbodies[0]["fileName"], 	filename); 
      NastranReader reader(filename.c_str());
      *m = TriMesh_d(reader,flipnormals);
      //mesh.push_back (new SPH::TriMesh(reader,flipnormals ));  
    }

    // Verifica que la construcción fue exitosa
    if (m == nullptr || m->nodecount == 0) {
        cout << "ERROR: Mesh construction failed!" << endl;
        return 0;
    }

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
    cout << "new mesh nodecount "<<m->allocated_meshes<<endl;
    dom_d->addMeshData(*m);
    delete m;
    cout << "Done."<<endl;
    
    cout << "MESH NODE COUNT "<<msh->nodecount<<endl;
  }//added another rigid body mesh
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
  #ifdef CUDA_BUILD
  calcMinEdgeLengthKernel<<<1,1>>>(dom_d); //TODO: PARALLELIZE
  cudaDeviceSynchronize();
  double *d_value;
  cudaMalloc((void**)&d_value, sizeof(double));

  getMinLengthKernel<<<1,1>>>(dom_d, d_value);
  cudaDeviceSynchronize(); 
  cudaMemcpy(&dx, d_value, sizeof(double), cudaMemcpyDeviceToHost);
  #else
  dom_d->calcMinEdgeLength();
  dx = dom_d->getMinLength();
  double size = 0.0;
  if(!readValue(mesh["elemSize"], size)){
    cout << "Setting remesh size at MIN ELEM SIZE:"<<dx<<endl;
  dom_d->setRemeshLength(dx);
  } else {
    cout << "Setting remesh size MANUALLY AT "<< size<<endl;    
    dom_d->setRemeshLength(size);
  }
  #endif
  

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
  

  
  int fixcount =0;
  int velcount =0;
  int xyzfixcount[] = {0,0,0};

  for (int i=0;i<dom_d->getNodeCount();i++){
     
     //~ #ifdef CUDA_BUILD
     //~ #else
     //~ if (dom_d->getPosVec3(i).z <0.0005) {
        //~ for (int d=0;d<3;d++)
          //~ dom_d->AddBCVelNode(i,d,0);
        //~ //dom_d->AddBCVelNode(i,2,0);
        //~ fixcount++;
        //~ //cout << "node "<< i<<" fixed "<<endl;
      //~ }
     //~ #endif

    //~ #ifdef CUDA_BUILD
    //~ #else
    //~ if (abs(dom_d->getPosVec3(i).x)<0.001 && abs(dom_d->getPosVec3(i).y<0.001 ) ) {
       //~ for (int d=0;d<2;d++)
         //~ dom_d->AddBCVelNode(i,d,0);
       //~ //dom_d->AddBCVelNode(i,2,0);
       //~ fixcount++;
       //~ //cout << "node "<< i<<" fixed "<<endl;
     //~ }
    //~ #endif
    
    //~ #ifdef CUDA_BUILD
    //~ #else
    //~ if (dom_d->getPosVec3(i).x <0.0004) {
       //~ dom_d->AddBCVelNode(i,0,0);
       //~ fixcount++;
       //~ //cout << "node "<< i<<" fixed "<<endl;
     //~ }
    //~ #endif
    
        //~ #ifdef CUDA_BUILD
    //~ #else
    //~ if (dom_d->getPosVec3(i).y <0.0004) {
       //~ dom_d->AddBCVelNode(i,1,0);
       //~ fixcount++;
       //~ //cout << "node "<< i<<" fixed "<<endl;
     //~ }
    //~ #endif
    
    //~ #ifdef CUDA_BUILD
    //~ #else    
     //~ if (dom_d->getPosVec3_h(i).z > 0.03-0.0005 ) {
       //~ dom_d->AddBCVelNode(i,0,-0.0);
       //~ dom_d->AddBCVelNode(i,1,-0.0);
       //~ dom_d->AddBCVelNode(i,2,-1.2);
       //~ cout << "Node "<<i <<" vel "<<endl;
       //~ velcount++;
     //~ }     
    //~ #endif
    
    for (int d=0;d<3;d++){
      if (xyzsym[d]){
        #ifdef CUDA_BUILD
        #else
          double coord;
          if      (d==0)  coord = dom_d->getPosVec3(i).x;
          else if (d==1)  coord = dom_d->getPosVec3(i).y;
          else            coord = dom_d->getPosVec3(i).z;
          if (coord < symtol ) {
          dom_d->AddBCVelNode(i,d,0);
          xyzfixcount[d]++;
                
          }
        #endif
      }
    }

    
  }
  //initElemArrayCPU (this,sigma_y,1,300.0e6)  


  
  cout << "FIXED "<<fixcount<< " NODES"<<endl;  
  cout << "VEL  "<<velcount<< " NODES"<<endl;  
  
  //AFTER THIS CALL
  dom_d->AllocateBCs();
  
 
  //
  cout << "Reading contact "<<endl;
  #ifdef CUDA_BUILD
  //cudaMalloc((void**)&dom_d->trimesh,         rigbodies.size()* sizeof(SPH::TriMesh_d));
  #else
    
  #endif
    
	cout << "Element Count "<<dom_d->getElemCount()<<endl;

    
    
  //dom_d->out_file.close(); // If done writing    
    

  }//if Json

	dom_d->SolveChungHulbert ();
  
  
	cout << "Program ended."<<endl;
      
	} else { //ARG >1
    cout << "Please provide an input file."<<endl;
  }
  
	
}
