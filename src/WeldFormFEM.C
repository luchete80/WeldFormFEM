#include "Domain_d.h"

#include <iostream>
#include "defs.h"

#include "lsdynaReader.h"
#include "nlohmann/json.hpp"
#include "Mesh.h"
#include <fstream>
#include <iomanip>	//ONY FOR GCC!!
#include "Input.h"

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


#endif

using namespace LS_Dyna;

int main(int argc, char **argv) {
  int dim = 3;
	if (argc > 1){
		string inputFileName=argv[1];	
		//std::ifstream i(argv[1]);    
  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);

	Domain_d *dom_d;

  #ifdef CUDA_BUILD
	report_gpu_mem();
	gpuErrchk(cudaMallocManaged(&dom_d, sizeof(MetFEM::Domain_d)) );
	report_gpu_mem();
  #else
    dom_d = new Domain_d;
  #endif    
    
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
    readValue(config["outTime"], out_time);
    readValue(config["simTime"], sim_time);

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
  E   = 70.0e9;
  nu  = 0.3;
  rho = 2700.0;
  

  cout << "Setting density"<<endl;
  dom_d->setDensity(rho);
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
    // cout << "Material Constants, Et: "<<c[0]<<endl;
    // material_h->Material_model = BILINEAR;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Bilinear )); //
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Bilinear), cudaMemcpyHostToDevice);	

    dom_d->AssignMaterial(material_h);
  } 
  cout << "Done."<<endl;
  // else if (mattype == "Hollomon")    {
    // // material_h  = new Hollomon(el,Fy,c[0],c[1]);
    // // cout << "Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
    // // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Hollomon));
    
    // material_h  = new Material_(el);
    // material_h->InitHollomon(el,Fy,c[0],c[1]);
    // material_h->Material_model = HOLLOMON;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Material_));
    
    // //init_hollomon_mat_kernel<<<1,1>>>(dom_d); //CRASH
    // //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Hollomon*), cudaMemcpyHostToDevice);	
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Material_), cudaMemcpyHostToDevice);	 //OR sizeof(Hollomon)??? i.e. derived class
    
  
  // } else if (mattype == "JohnsonCook") {
    // //Order is 
                               // //A(sy0) ,B,  ,C,   m   ,n   ,eps_0,T_m, T_transition
   // //Material_ *material_h  = new JohnsonCook(el,Fy, c[0],c[1],c[3],c[2],c[6], c[4],c[5]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
    
    // //Only 1 material to begin with
    // //cudaMalloc((void**)&dom_d->materials, 1 * sizeof(JohnsonCook ));
    // //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(JohnsonCook), cudaMemcpyHostToDevice);	
    // cout << "Material Constants, B: "<<c[0]<<", C: "<<c[1]<<", n: "<<c[2]<<", m: "<<c[3]<<", T_m: "<<c[4]<<", T_t: "<<c[5]<<", eps_0: "<<c[6]<<endl;
  // } else                              printf("ERROR: Invalid material type.


      
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

  if (contact){
    cout << "Searching external nodes"<<endl;
    dom_d->SearchExtNodes();
    
    TriMesh_d *msh = new TriMesh_d();
    
    //AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
    cout <<"Creating plane mesh..."<<endl;
   //void TriMesh_d::AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens)
    msh->AxisPlaneMesh(2, false, make_double3(0.095,-0.005,0.024), make_double3 (.115,0.0015,0.024),  2);
    msh->CalcSpheres();  //NFAR Done Once if mesh is rigid
    cout <<"Done"<<endl;
    msh->SetVel(make_double3(0.0,0.,-0.48));
    dom_d->setTriMesh(msh);
    
    dom_d->setContactOn();
  }

	double dt = 0.7 * dx/(mat_cs);
  //double dt = 0.800e-5;
  dom_d->SetDT(dt); 
  cout << "End Time: "<<sim_time<<endl;
  dom_d->SetEndTime (sim_time);
  dom_d->setdtOut(out_time);
  //dom_d->SetEndTime (1000.0*dt);
  

  int fixcount =0;
  int velcount =0;
  for (int i=0;i<dom_d->getNodeCount();i++){
    
    if (dom_d->getPosVec3(i).z <0.025) {
      for (int d=0;d<3;d++)dom_d->AddBCVelNode(i,d,0);
      fixcount++;
      cout << "node "<< i<<" fixed "<<endl;
    }
    
    if (dom_d->getPosVec3(i).z > 0.616-0.025 ) {
      dom_d->AddBCVelNode(i,0,-0.0);
      dom_d->AddBCVelNode(i,1,-0.0);
      dom_d->AddBCVelNode(i,2,-1.0);
      cout << "Node "<<i <<" vel "<<endl;
      velcount++;
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
  cudaMalloc((void**)&dom_d->trimesh,         rigbodies.size()* sizeof(SPH::TriMesh_d));
  #else
    
  #endif
    
	cout << "Element Count "<<dom_d->getElemCount()<<endl;

    
    
    
    

  }

	dom_d->SolveChungHulbert ();
  
  
	cout << "Program ended."<<endl;
      
	} else {
    cout << "Please provide an input file."<<endl;
  }
	
}
