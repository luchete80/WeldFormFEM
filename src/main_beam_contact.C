#include "Domain_d.h"

#include <iostream>
#include "defs.h"

#include "Mesh.h"

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

int main(){

  int dim = 3;
  double tf = 1.0e-3;
  
	Domain_d *dom_d;

  #ifdef CUDA_BUILD
	report_gpu_mem();
	gpuErrchk(cudaMallocManaged(&dom_d, sizeof(MetFEM::Domain_d)) );
	report_gpu_mem();
  #else
    dom_d = new Domain_d;
  #endif
	
  double3 V = make_double3(0.0,0.0,0.0);
  double dx = 0.006;
  double Lx = 0.1	;
  double Ly = 0.024;
  double Lz = 0.012;
	double3 L = make_double3(Lx,Ly,Lz);
  
  //dom_d->Nproc = 8;
  dom_d->setNproc(8);
  
	double r = dx/2.0;
	
	dom_d->AddBoxLength(V,L,r,true);
 

  ////// MATERIAL  
  double E, nu, rho;
  E   = 70.0e9;
  nu  = 0.3;
  rho = 2700.0;
  

  
  dom_d->setDensity(rho);
  
  cout << "Creating Material..:"<<endl;
  Material_ *mat_h = (Material_ *)malloc(dom_d->getElemCount() * sizeof(Material_ *)); 
  Elastic_ el(E,nu);
  // cout << "Mat type  "<<mattype<<endl;

  Material_ *material_h;
  double Ep, c[6];
  // MATERIAL
  //TODO: MATERIALS SHOULD BE A VECTOR

  
  // !ONLY FOR ELASTIC
  // dom%mat_E = young
  // dom%mat_nu = poisson
  
  double mat_modK= E / ( 3.0*(1.0 -2.0*nu) );
  double mat_G= E / (2.0* (1.0 + nu));
  
  // dom%mat_K = mat_modK !!!TODO CREATE MATERIAL
  
  double mat_cs = sqrt(mat_modK/rho);
  // mat_cs0 = mat_cs
  // print *, "Material Cs: ", mat_cs
  
  // elem%cs(:) = mat_cs
  
  // dt = 0.7 * dx/(mat_cs)
  
  string mattype = "Bilinear";
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

	double dt = 0.7* dx/(mat_cs);
  //double dt = 0.800e-5;
  dom_d->SetDT(dt); 
  dom_d->SetEndTime (1.e-3);
  //dom_d->SetEndTime (1000.0*dt);
  
  
  TriMesh_d *msh = new TriMesh_d();
  
  //AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
  cout <<"Creating plane mesh..."<<endl;
  msh->AxisPlaneMesh(2, true, make_double3(0,0,0), make_double3 (1,1,0),  2);
  cout <<"Done"<<endl;
  //msh->SetVel(make_double3(0.0,0.,-1.0));
  //dom_d->setTriMesh(msh);
  
  //dom_d->SearchExtNodes();
  
  int fixcount =0;
  int velcount =0;
  for (int i=0;i<dom_d->getNodeCount();i++){
    
    if (dom_d->getPosVec3(i).x <r/2.0) {
      for (int d=0;d<3;d++)dom_d->AddBCVelNode(i,d,0);
      fixcount++;
      cout << "node "<< i<<" fixed "<<endl;
    }
    
    if (dom_d->getPosVec3(i).x > (Lx - 1.5*r) && dom_d->getPosVec3(i).y > (Ly -1.5*r) ) {
      dom_d->AddBCVelNode(i,1,-0.48);
      cout << "Node "<<i <<" vel "<<endl;
      velcount++;
    }      
    
  }
  cout << "FIXED "<<fixcount<< " NODES"<<endl;  
  cout << "VEL  "<<velcount<< " NODES"<<endl;  
  
  //AFTER THIS CALL
  dom_d->AllocateBCs();
  
 

	cout << "Element Count "<<dom_d->getElemCount()<<endl;
	dom_d->SolveChungHulbert ();
  
  
	cout << "Program ended."<<endl;
  
  //dom_d->WriteToCSV("test.csv");
	
	
}
