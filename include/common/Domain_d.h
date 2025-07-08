#ifndef _DOMAIN_D_CUH_
#define _DOMAIN_D_CUH_

//////////////////// COMMON DOMAIN FOR CPU/GPU ///////////////////

#include "defs.h"

#include <stdio.h>

#ifdef  CUDA_BUILD
#include <cuda.h>
#endif

#include "utils.h"

#include "Material.cuh"

#include <vector>
#include <string>
#include "lsdynaReader.h"

#include "parallel_for_each.h"

#include <fstream>

#define _QUA2D_ 0
#define _TRI2D_ 1
#define _TET2D_ 2
     
//Assuming calling from Domain_d *dom   
//You can add dom as another macro argument

//FROM DOMAIN
#define initNodalArrayCPU(dom,a,dim,val) for (int n=0;n<dom->m_node_count;n++){\
                                        for (int d=0;d<dim;d++)\
                                        a[n*dim + d] = val;\
                                      }

#define initElemArrayCPU(dom,a,dim,val) for (int n=0;n<dom->m_elem_count;n++){\
                                        for (int d=0;d<dim;d++)\
                                        a[n*dim + d] = val;\
                                      }

#if CUDA_BUILD

#define initElemArray (dom,a,dim,val) for (int n=0;n<dom->m_elem_count;n++){\
                                        for (int d=0;d<dim;d++)\
                                        a[n*dim + d] = val;\
                                      }
#endif
                                      
//// THIS IS IF THE VARIABLE IS IN THE GAUSS POINTS (ASSUMING 1 GAUSS POINT)  
//// IS NOT THE SAME FOR ELEMENT NODAL VARIABLE (LIKE ELEMENT FORCES)        
////FOR ELEMENT FORCES IT WOULD BE:  ///a[n*dim + d]+=v[offset + ne*dim + d];\
                          
                              //FOR INIT CALL PREVIOUS ARRAY
#define avgScalar(v,a,dim)    for (int n=0;n<dom->m_node_count;n++){\
                                for (int d=0;d<dim;d++)\
                                  a[n*dim + d] = 0.0;\
                                for (int e=0; e<dom->m_nodel_count[n];e++) {\
                                  int eglob   = dom->m_nodel     [dom->m_nodel_offset[n]+e];\
                                  int ne      = dom->m_nodel_loc [dom->m_nodel_offset[n]+e];\
                                  int offset  = eglob * dom->m_nodxelem * dim;\
                                  for (int d=0;d<dim;d++) a[n*dim + d]+=v[eglob*dim + d];\
                                  }\
                                for (int d=0;d<dim;d++) a[n*dim + d]/=dom->m_nodel_count[n];\
                              }

enum dom_type {_Plane_Strain_=0,_Plane_Stress_, _Axi_Symm_, _3D_};


class Matrix;

class BC_Node {
  public:
  BC_Node(const int &nod, const int &dof, const int &val):
  m_node(nod){
    
  }
  int     m_node;
  int     m_dof;
  double  m_val;
};

namespace MetFEM{
  
enum DomainType { IMPL, EXPL };
 

// Structure to define a face with 4 nodes
//FOR HEXA IS 4
//PARALLELIZE WITH GPU : TODO
//// ORIGINALLY MEANS HEXA
#define ELNOD  4   //ORIGINALLY 8
#define FACENOD 3  //ORIGINALLY 4
#define ELFAC  4   //ORIGINALLY 6
//OLD FOT HEXA, CHANGE IT

struct Face {
    int nodes[FACENOD];
    int count; // Number of occurrences of this face
    int elem_id;  // <- store the originating element
    int other_elem;   // second elem

};


enum BC_TYPE {Velocity_BC=0, Force_BC, Temperature_BC, Convection_BC, Symmetry_BC};

struct boundaryCondition {
	int 	zoneId;
	BC_TYPE 	type;	// ENUM TYPE Velocity, Force, Temperature
	bool 	free;	//is necessary??
	bool  is_type[4];  ////ACCORDING TO BC_TYPE {Velocity_BC=0, Force_BC, Temperature_BC, Convection_BC};
	int 	valueType;		//0: Constant, 1 amplitude table
	double3 value;       //If constant
	double3 value_ang;       //Angular value
	int 	ampId;			//if valuetype == 1
	double 	ampFactor;		//if valuetype == 1
  double cv_coeff;
  double T_inf;
  double T;
};

class TriMesh_d;
class VTKWriter;
class VTUWriter;
class ReMesher;
class ReMesher;
class Solver;
class Solver_Eigen;

class Domain_d {
  friend class VTKWriter;
  friend class VTUWriter;
  friend class ReMesher;
  friend class Solver;
  friend class Solver_Eigen;
public:
  Domain_d (std::string);
  Domain_d (){
    m_axisymm_vol_weight = false;
    m_domtype = _Plane_Strain_;
    bc_count[0]=bc_count[1]=bc_count[2]=0;
    contact = false;
    m_thermal = false;
    m_remesh_interval = 1e10;
  }
  void setNproc(const int &n){Nproc=n;}
  void SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  void SetDimensionImplicit(const int &node_count, const int &elem_count); //TO NOT USE virtual (GPU issues)
  void AddBoxLength(vector_t const & V, vector_t const & L, const double &r, const bool &red_int = true, const bool &tritetra = false);
  void CreateFromLSDyna(LS_Dyna::lsdynaReader &reader);
  ///// (CUDA HOST) FUNCTIONS 
  inline dev_t double2 getPosVec2(const int &n){
    return make_double2(x[m_dim*n], x[m_dim*n+1]);
    };
  inline dev_t vector_t getPosVec3  (const int &n){return make_vector_t(x[m_dim*n], x[m_dim*n+1], x[m_dim*n+2]);}; //the same
  

  #ifdef CUDA_BUILD
  inline vector_t getPosVec (int &n){return make_vector_t(x_h[m_dim*n], x_h[m_dim*n+1], x_h[m_dim*n+2]);};
  inline double2  getPosVec2_h(const int &n){return make_double2 (x_h[m_dim*n], x_h[m_dim*n+1]);}; //the same
  inline vector_t getPosVec3_h(const int &n){return make_vector_t(x_h[m_dim*n], x_h[m_dim*n+1], x_h[m_dim*n+2]);}; //the same
  #else
  inline dev_t vector_t getPosVec3_h(const int &n){return make_vector_t(x[m_dim*n], x[m_dim*n+1], x[m_dim*n+2]);}; //the same of getPosVec3
  inline dev_t double2 getPosVec2_h(const int &n){
    return make_double2(x[m_dim*n], x[m_dim*n+1]);
    };
  #endif

  void setTriMesh(TriMesh_d *m){trimesh = m;}
  void addMeshData(const TriMesh_d &m);
  
  
  dev_t void SearchExtNodes();

  inline vector_t getAccVec (const int &n){return make_double3(a[m_dim*n], a[m_dim*n+1], a[m_dim*n+2]);};  
  inline vector_t getContForceVec(const int &n){return make_vector_t(contforce[m_dim*n], contforce[m_dim*n+1], contforce[m_dim*n+2]);}
  
  #ifdef CUDA_BUILD
  inline vector_t getDispVec(const int &n){return make_vector_t(u_h[m_dim*n], u_h[m_dim*n+1], u_h[m_dim*n+2]);};
  inline vector_t getVelVec (const int &n){return make_vector_t(v_h[m_dim*n], v_h[m_dim*n+1], v_h[m_dim*n+2]);};
  inline vector_t getIntForceVec(const int &n){/*return make_vector_t(m_fi[m_dim*n], m_fi[m_dim*n+1], m_fi[m_dim*n+2]);*/};
  #else
  //inline vector_t getAccVec (const int &n){return make_vector_t(a[m_dim*n], a[m_dim*n+1], a[m_dim*n+2]);};
  inline vector_t getVelVec (const int &n){return make_vector_t(v[m_dim*n], v[m_dim*n+1], v[m_dim*n+2]);};
  inline vector_t getDispVec(const int &n){return make_vector_t(u[m_dim*n], u[m_dim*n+1], u[m_dim*n+2]);};
  inline vector_t getIntForceVec(const int &n){return make_vector_t(m_fi[m_dim*n], m_fi[m_dim*n+1], m_fi[m_dim*n+2]);};
  //inline vector_t getContForceVec(const int &n){return make_vector_t(contforce[m_dim*n], contforce[m_dim*n+1], contforce[m_dim*n+2]);}
  #endif
  
  dev_t void printVec(double*);
  dev_t void printSymmTens( double *);
  
  void setNodElem(int *elnod);//For assembly and parallel processing
  
  void WriteToVTK(char *);
  int WriteToCSV(char *);
  
  dev_t void calcElemJacobian ();
  dev_t void calcElemJAndDerivatives/*_FullInt*/();
  dev_t void calcElemMassMat();

	dev_t void calcElemStrainRates();
  dev_t void calcElemForces();
  dev_t void calcElemHourglassForces();
  dev_t void calcElemPressure(); //FROM STRAIN
  dev_t void smoothPressureField(double gamma); //Laplace Smooth
  dev_t void calcElemPressureFromJ();
  dev_t void calcElemPressureANP(); //AVERAGE NODAL POINT
  dev_t void calcElemPressureElementBased();
  dev_t void calcElemPressureANP_Nodal();
  dev_t void calcElemPressureANP_Nodal_HG();
  
  dev_t void calcElemPressure_Hybrid();
  dev_t void calcElemPressure_Hybrid_VolHG();

  dev_t void assemblyMassMatrix();
  dev_t void assemblyForces();
  
  host_   void AddBCVelNode(const int &node, const int &dim, const double &val);
  host_   void AllocateBCs();
  
  void calcMassDiagFromElementNodes(const double &rho); // To use existing array, NOT IN PARALLEL

  dev_t void ImposeBCA(const int dim); /// DO NOT USE REFERENCESSS!!!!!!
  host_ void ImposeBCAAllDim();

  dev_t void ImposeBCV(const int dim); /// DO NOT USE REFERENCESSS!!!!!!
  host_ void ImposeBCVAllDim();
  
  dev_t void calcMinEdgeLength();
  
  ///// ATENTION! THIS IS Deriv x DETJ
  inline dev_t double & getDerivative(const int &e, const int &gp, const int &i, const int &j); //I AND J ARE: DIMENSION AND NODE
  inline dev_t void     setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &); //I AND J ARE: DIMENSION AND NODE
  inline dev_t void     setDerivative(const int &e, const int &gp, Matrix *); //I AND J ARE: DIMENSION AND NODE
  
  void setDensity(const double &r);
  
  // template <typename T>
  // void setVector(,const double &r);
  
	int threadsPerBlock, blocksPerGrid; //TO BE USED BY SOLVER
	
	const int & getElemCount()const{return m_elem_count;}
  unsigned int & getElemNode(const int &e, const int &n)const{return m_elnod[m_nodxelem*e+n];}
	const int & getNodeCount()const{return m_node_count;}
  vector_t getNodePos3(const int &n){
    #ifdef CUDA_BUILD
    return make_vector_t (x_h[m_dim*n],x_h[m_dim*n+1],x_h[m_dim*n+2]);
    #else
    return make_vector_t (x  [m_dim*n],x  [m_dim*n+1],x  [m_dim*n+2]);      
    #endif
  }


  double2 getNodePos2(const int &n){
    #ifdef CUDA_BUILD
    return make_double2 (x_h[m_dim*n],x_h[m_dim*n+1]);
    #else
    return make_double2 (x  [m_dim*n],x  [m_dim*n+1]);      
    #endif
  }

  
  void AssignMaterial (Material_ *material_h); //Create and copy material
  
  dev_t void AssignMatAddress(); //Assign particle data to material array to zero array
  void AssignMatAddress_(); //Assign particle data to material array to zero array
  
  dev_t void CalcStressStrain(double dt);
  dev_t void Calc_Elastic_Stress (const double dt);
  
  ///// AXISYMM VARS
  dev_t void Calc_Element_Radius(); //For axisymm
  dom_type m_domtype;
  bool m_axisymm_vol_weight;
  void setDomtype();

  //TODO: CHANGE THIS
  inline dev_t double & getSigma  (const int e, const int gp, int i, int j){int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}};/*printf("i j symm idx %d %d %d\n",i,j,symm_idx[i][j]);*/if (j<3) return m_sigma   [e*m_gp_count*6 + symm_idx[i][j]];}
  inline dev_t double & getStrRate(const int e, const int gp, int i, int j){int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}}; if (j<3) return m_str_rate[e*m_gp_count*6 + symm_idx[i][j]];}
  
  dev_t void CalcElemInitialVol();
  
  dev_t void InitUVA(); 
  
  dev_t void CalcElemVol();
  //Specially for tetra
  dev_t void CalcNodalVol();
  dev_t void CalcNodalMassFromVol();
  
  dev_t void calcElemDensity();
    
  dev_t void calcAccel();
  dev_t void UpdatePrediction();
  dev_t void UpdateCorrectionAccVel();
  dev_t void UpdateCorrectionPos();
  dev_t void calcTotMass();
  dev_t void InitElemValues(double *arr, double val = 0.0);
  dev_t void InitStressesFromMat();
    
  inline dev_t double getRadius(int e, int gp){return m_radius[e*m_gp_count+gp];}
	
  //__device__ vector_t & getVElem(const int &e, const int &n){return v[m_elnod[e*m_nodxelem+n]];}
  inline dev_t double  getVElem(const int &e, const int &n,const int &d){return v[m_dim*m_elnod[m_nodxelem*e+n]+d];}  
  
  inline dev_t vector_t getV(const int &n){return make_vector_t(v[m_dim*n], v[m_dim*n+1], v[m_dim*n+2]);}  
  
  
	void SolveChungHulbert();
  void setdtOut(const double &t){m_dtout=t;}
  void SetDT(const double &dt_){dt=dt_;}
  void SetEndTime(const double &tf_){end_t=tf_;}
	int & getDim(){return m_dim;}
  void setProcCount(const int &np){Nproc = np;} //for CPU only
  
    
  void dev_t CalcContactForces(); //standard, not predicted
  void setContactOn(){contact = true;}
  bool isContactOn(){return contact;}
  TriMesh_d* getTriMesh(){return trimesh;}
  void InitValues();
  double dev_t getMinLength(){return m_min_length;}
  
  dev_t void calcNodalPressureFromElemental();
  
  bool m_auto_contact; // if not define external nodes
  void setRemeshInterval(int i) {m_remesh_interval = i;}
  void setRemeshLength(const double &d){m_remesh_length=d;}
  void Free();
  void calcArtificialViscosity();
  
  /////// THERMAL
  void ThermalCalcs();
  void setThermalOn(){m_thermal = true;}
  void calcInelasticHeatFraction();
  void setTemp(const double &val){
    #ifdef  CUDA_BUILD
    
    #else
    for (int i=0;i<m_node_count;i++)T[i]=val;
    #endif
    
    }
  void calcContactForceFromPressure();
  void setNode(const int &i, const double &_x, const double &_y, const double &_z);
    
  //--------------------------------------------------------------------------------------------------------------------------------
  
  void CalcExtFaceAreas();
  std::ofstream out_file;
  void setCFL(const double &f){m_cfl_factor = f;}
  
    // IMPLICIT FUNCTIONS
  //--------------------------------------------------------------------------------------------------------------------------------
  Matrix getElemBMatrix(const int &e);
  void CalcMaterialStiffElementMatrix();
  void CalcGeomStiffElementMatrix();  
  void CalcElemIntForces();
	void Solve();
	void ElasticSolve();
   
  //--------------------------------------------------------------------------------------------------------------------------------

  
  void setFixedDt(const bool &f){m_fixed_dt = f;}
  const bool &ifFixedDt()const{return m_fixed_dt;}
  double m_artifvisc[2];
  
protected:
  double          m_tot_mass; //Only for testing, not needed
  dev_t int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
	int 						m_dim;
  int             m_node_count, m_elem_count;
	
	unsigned int 		*m_elnod;              /// TO SHIFT 
  
  
  unsigned int   *m_elnodoffset; //FROM FIRST ELEMENT NODE
	
	double 				 *x, *x_h; //Vector is double, x_h is used only on host
	double 				 *v; //CHANGED TO DOUBLE
	double 				 *a, *prev_a;
	double         *u, *u_dt;
  double         *v_h,*u_h, *a_h;
  double         *m_mglob, *m_mdiag, *m_ematm; //Global Matrix (TODO: MAKE SPARSE) and Diagonal masses, and element masses
  
	double 					*p, *rho, *rho_0;
  
  double          *ut_prev;   //FOR HISTORICAL CONTACT
  
  //////Thermal
  bool m_thermal;
  double *T;   /// temperature
  double *m_dTedt; //elem node
  double *ps_energy;
  double *q_cont_conv;
  double *node_area; //THERMAL CONTACT
  double  contact_hc; //MOVETO CONTACT INTERFACE
  double *m_elem_length; //FOR CONTACT
  double *p_node; //FOR CONTACT AND ANP (ANP will be overriden by HG later)
  double *m_elem_area;
  
  int 	 *m_mesh_in_contact; //ID (Pos) with mesh in contact
  
  double m_min_length;
  double m_min_height;
  double m_cfl_factor;
  int             m_remesh_interval;
  double          m_remesh_length;

  double          *pl_strain, *sigma_y;
  
  
  double          dt, end_t;
  
  double          m_alpha, m_beta, m_gamma;
  double          *vol_0, *vol;  //Element initial and current volume
  
  double          *m_voln,*m_vol_0n; //Nodal values, useful for TETRA nodal Average Nodal Pressure
  
  //TO NOT OCCUPY SO MUCH MEMORY
  int             *bcx_nod, *bcy_nod, *bcz_nod;
  double          *bcx_val, *bcy_val, *bcz_val;
  int             bc_count[3];
  
  ////// PARALLELIZATION
  int             Nproc;
  
  std::vector<double>     bcx_val_h,bcy_val_h,bcz_val_h;
  std::vector<int>        bcx_nod_h,bcy_nod_h,bcz_nod_h;
  
  Material_ **mat; //pointer to material of each particle
  Material_ *materials; //WAll materials 
  
  bool            m_red_int;  //Reduced integration, 1 GAUSS POINT
  int             m_gp_count; //Assuming constant gauss points
  int             m_nodxelem; //THIS IS INTENDED TO BE MODIFIED BY m_nodxelem_e which is a matrix
  int             *m_nodxelem_e;
  //PLASTIC HEAT 
  double          *m_strain_pl_incr; //Increment of plastic strain, useful for plastic heat generation.
  double          *m_q_plheat;
  double_t        m_pl_energy;
  
  /////// ELEMENT VARIABLES
  /////////////////////// //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION ///////////////////////////////////
  Matrix          **m_dHrs;     //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  Matrix          **x2;         //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  Matrix          **dH_dxyz; 
  
  Matrix          *m_jacob;
  int             *elnod_h; ////// USED TO COMPUTE GLOBAL M MATRIX WHICH IS COMPUTED IN CPU (TODO: CHANGE)       
  double          *m_str_rate, *m_rot_rate;
  double          *m_f_elem, *m_f_elem_hg;    //ELEMENT
  double          *m_fi, *m_fe; //NODAL
  double          *m_sigma, *m_tau;
	double          *m_radius;
  //Updated lagrangian formulation
  //real(fp_kind), dimension(:,:,:,:), allocatable :: BL,BNL, jacob, dHxy,dHxy_detJ, dHxy0,math, dHrs !!!DIM: e,gp,,:,:, is it necesary to store dHrs??? is only because is used twice, at J and dHxy
  double         *dHxy_detJ ; ////NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
  double         *m_detJ;
  
  ////// THESE ARE SUBDIVIDED FOR HIGHER ACCESS SPEED (LOWER OFFSET)
  double         *m_dH_detJ_dx ; //GAUSS_POINT, AGRUPATED BY ELEM 1 GP1, ELEM 1 GP 2 .... ELEM 2 GP 1..
  double         *m_dH_detJ_dy ; 	
  double         *m_dH_detJ_dz ; 

  double         *m_H;          // Shape functions, for mass calculation 
	
  double					Time;    				//Current time of simulation at each solving step
	double					dtmin;			//Minimum Time Step
	double					dtint;			//Initial Time Step
  double          m_dtout;  //VTK outtime  
  
  ////////////////////////////// REDUCTION THINGS ////////////////////////////////////////
  //THESE ARE USED IF NODAL FORCES WANTED TO BE REDUCED (SUMMED) NOT BY ATOMIC (GPU) OR LOCK (OPENMP CPU)
  bool            m_nonlock_sum;
  int            *m_nodel;           // NODE ELEMENTS: ELEMENTS SHARED BY EACH NODE [nodxelem* node_count] call:  m_nodel[nod,elcount] =EL i,e m_nodel[nod_offset+elcount]
  int            *m_nodel_loc;        //
  int            *m_nodel_offset;    //OFFSET OF THE
  int            *m_nodel_count;    

  //////////////////////////// USED FOR CONTACT
  unsigned int    *m_elnod_count;   /// FOR CONTACT, TO REPLACE FOR m_node_count
  unsigned int    *m_contsurf_count;
  unsigned int    *m_contsurf_elemcount;   //FOR EACH OF THE ABOVE  
  unsigned int    *m_contsurf_elem;        //ELEMENT POS OF THE CONTACT ELEMENT 

  Face *faceList;
  int m_faceCount;
  
  int *m_elem_neigh;
  int *m_elem_neigh_count;
    
  ////////////////////// CONTACT 
	// TODO, EACH RIGID PARTICLE SHOULD 
  int   *contelem; //ELEMENT OF TRIMESH FROM "RIGID" PARTICLE, ALL FIRST PARTICLES ARE ZERO
  TriMesh_d *trimesh;
  int trimesh_count;
  int *mesh_id; //particle mesh ID	
	bool *ext_nodes;
  int ext_nodes_count;
  double *contforce; 
  bool contact;
  
  bool m_fixed_dt;
  
  /////// IMPLICIT
  //////////////////////////////// IMPLICIT THINGS
  Matrix **m_Kmat;   //MATERIAL PART OF 
  Matrix **m_Kgeo;   //MATERIAL PART OF   
  
  //// DEFORMATION GRADIENT
  Matrix **m_Fel;    //def
  
  
  Solver* m_solver;

};

#ifdef  CUDA_BUILD

__global__ void ImposeBCVKernel(Domain_d *dom_d, int d);
__global__ void ImposeBCAKernel(Domain_d *dom_d, int d);

__global__ void calcElemJAndDerivKernel(Domain_d *dom_d);
__global__ void calcElemMassMatKernel(Domain_d *dom_d);
__global__ void calcElemStrainRatesKernel(Domain_d *dom_d);

__global__ void assemblyForcesKernel(Domain_d *dom_d);
__global__ void assemblyMassMatrixKernel(Domain_d *dom_d);

__global__ void calcElemForcesKernel          (Domain_d *);
__global__ void calcElemHourglassForcesKernel (Domain_d *);

__global__ void calcElemPressureKernel (Domain_d *);
__global__ void calcElemDensityKernel  (Domain_d *);

__global__ void UpdatePredictionKernel(Domain_d *dom_d);
__global__ void UpdateCorrectionAccVelKernel(Domain_d *dom_d);
__global__ void UpdateCorrectionPosKernel   (Domain_d *dom_d);

__global__ void calcStressStrainKernel(Domain_d *dom_d, double dt);

__global__ void calcElemVolKernel(Domain_d *dom_d);
//__global__ void calcNoadlVolKernel(Domain_d *dom_d);

__global__ void calcElemInitialVolKernel(Domain_d *dom_d);

__global__ void AssignMatAddressKernel(Domain_d *dom/*, Material_ *mat*/);

__global__ void calcTotMassKernel(Domain_d *dom_d);

__global__ void calcMassDiag(Domain_d *dom_d);

__global__ void calcAccelKernel(Domain_d *dom_d);

__global__ void printVecKernel(Domain_d *dom_d, double *);
__global__ void printSymmTensKernel(Domain_d *dom_d,  double *);
__global__ void calcMinEdgeLength(Domain_d *dom_d);

__global__ void InitElemValuesKernel(Domain_d *dom_d, double *arr, double val = 0.0);
__global__ void InitStressesFromMatKernel(Domain_d *dom_d);
__global__ void CalcNodalVolKernel        (Domain_d *dom_d);
__global__ void CalcNodalMassFromVolKernel(Domain_d *dom_d);
__global__ void calcMinEdgeLengthKernel(Domain_d *dom_d);
__global__ void getMinLengthKernel(Domain_d *dom_d, double *d);

// Kernel to retrieve only the private value
__global__ void getMinLengthKernel(Domain_d *dom_d,  double *d_value);
#endif

//i: dim, j: node
inline dev_t double & Domain_d::getDerivative(const int &e, const int &gp, const int &i, const int &j){ //I AND J ARE: DIMENSION AND NODE
      if (i == 0)     return m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_nodxelem + j];
      else if (i==1)  return m_dH_detJ_dy[e*(m_nodxelem * m_gp_count) + gp * m_nodxelem + j];
      else if (i==2)  return m_dH_detJ_dz[e*(m_nodxelem * m_gp_count) + gp * m_nodxelem + j];
      else printf ("ERROR: WRONG DERIVATIVE DIMENSION.");
}

inline dev_t void Domain_d::setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &v){ //I AND J ARE: DIMENSION AND NODE
      if (i == 0)     m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_nodxelem + j] = v;
      else if (i==1)  m_dH_detJ_dy[e*(m_nodxelem * m_gp_count) + gp * m_nodxelem + j] = v;
      else if (i==2)  m_dH_detJ_dz[e*(m_nodxelem * m_gp_count) + gp * m_nodxelem + j] = v;
      else printf ("ERROR: WRONG DERIVATIVE DIMENSION.");
}

inline dev_t void Domain_d::setDerivative(const int &e, const int &gp, Matrix *m){ //I AND J ARE: DIMENSION AND NODE
      // for (int j = 0;j<3;j++)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i] = v;
      // else if (i==1)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i] = v;
      // else if (i==2)  m_dH_detJ_dx[e*(m_nodxelem * m_gp_count) + gp * m_gp_count + i] = v;
      //else printf ("ERROR: WRONG DERIVATIVE DIMENSION.");
}



}; //Namespace

//#include "Mesh.h"
//#include "Contact.C"

#endif
