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
#include "../lib/LSDynaReader/src/lsdynaReader.h"

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

class VTKWriter;
class Domain_d {
  friend class VTKWriter;
public:
  Domain_d (std::string);
  Domain_d (){}
  void SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  void AddBoxLength(vector_t const & V, vector_t const & L, const double &r, const bool &red_int = true);
  void CreateFromLSDyna(const LS_Dyna::lsdynaReader &reader);
  ///// (CUDA HOST) FUNCTIONS 
  #ifdef CUDA_BUILD
  inline vector_t getPosVec(const int &n){return make_vector_t(x_h[3*n], x_h[3*n+1], x_h[3*n+2]);};
  #else
  inline vector_t getPosVec(const int &n){return make_vector_t(x[3*n], x[3*n+1], x[3*n+2]);};
  inline vector_t getIntForceVec(const int &n){return make_vector_t(m_fi[3*n], m_fi[3*n+1], m_fi[3*n+2]);};
  #endif

  #ifdef CUDA_BUILD
  inline vector_t getDispVec(const int &n){return make_vector_t(u_h[3*n], u_h[3*n+1], u_h[3*n+2]);};
  inline vector_t getVelVec (const int &n){return make_vector_t(v_h[3*n], v_h[3*n+1], v_h[3*n+2]);};
  #else
  inline vector_t getVelVec (const int &n){return make_vector_t(v[3*n], v[3*n+1], v[3*n+2]);};
  inline vector_t getDispVec(const int &n){return make_vector_t(u[3*n], u[3*n+1], u[3*n+2]);};
  inline vector_t getIntForceVec(const int &n){return make_vector_t(m_fi[3*n], m_fi[3*n+1], m_fi[3*n+2]);};
  #endif
  
  dev_t void printVec(double*);
  dev_t void printSymmTens( double *);
  
  void WriteToVTK(char *);
  int WriteToCSV(char *);
  
  dev_t void calcElemJacobian ();
  dev_t void calcElemJAndDerivatives/*_FullInt*/();
  dev_t void calcElemMassMat();

	dev_t void calcElemStrainRates();
  dev_t void calcElemForces();
  dev_t void calcElemHourglassForces();
  dev_t void calcElemPressure(); //FROM STRAIN

  dev_t void assemblyMassMatrix();
  dev_t void assemblyForces();
  
  host_   void AddBCVelNode(const int &node, const int &dim, const double &val);
  host_   void AllocateBCs();
  
  void calcMassDiagFromElementNodes(const double &rho); // To use existing array, NOT IN PARALLEL

  dev_t void ImposeBCA(const int dim); /// DO NOT USE REFERENCESSS!!!!!!
  host_ void ImposeBCAAllDim();

  dev_t void ImposeBCV(const int dim); /// DO NOT USE REFERENCESSS!!!!!!
  host_ void ImposeBCVAllDim();
  
  
  ///// ATENTION! THIS IS Deriv x DETJ
  inline dev_t double & getDerivative(const int &e, const int &gp, const int &i, const int &j); //I AND J ARE: DIMENSION AND NODE
  inline dev_t void     setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &); //I AND J ARE: DIMENSION AND NODE
  inline dev_t void     setDerivative(const int &e, const int &gp, Matrix *); //I AND J ARE: DIMENSION AND NODE
  
  void setDensity(const double &r);
  
  // template <typename T>
  // void setVector(,const double &r);
  
	int threadsPerBlock, blocksPerGrid; //TO BE USED BY SOLVER
	
	const int & getElemCount()const{return m_elem_count;}
	const int & getNodeCount()const{return m_node_count;}
  vector_t getNodePos(const int &n){
    #ifdef CUDA_BUILD
    return make_vector_t (x_h[3*n],x_h[3*n+1],x_h[3*n+2]);
    #else
    return make_vector_t (x  [3*n],x  [3*n+1],x  [3*n+2]);      
    #endif
  }
  
  void AssignMaterial (Material_ *material_h); //Create and copy material
  dev_t void AssignMatAddress(); //Assign particle data to material array to zero array
  dev_t void CalcStressStrain(double dt);
  dev_t void Calc_Elastic_Stress (const double dt);

  //TODO: CHANGE THIS
  inline dev_t double & getSigma  (const int e, const int gp, int i, int j){int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}};/*printf("i j symm idx %d %d %d\n",i,j,symm_idx[i][j]);*/if (j<m_dim) return m_sigma   [e*m_gp_count*6 + symm_idx[i][j]];}
  inline dev_t double & getStrRate(const int e, const int gp, int i, int j){int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}}; if (j<3) return m_str_rate[e*m_gp_count*6 + symm_idx[i][j]];}
  
  dev_t void CalcElemInitialVol();
  
  dev_t void InitUVA(); 
  
  dev_t void CalcElemVol();
  dev_t void calcElemDensity();
    
  dev_t void calcAccel();
  dev_t void UpdatePrediction();
  dev_t void UpdateCorrectionAccVel();
  dev_t void UpdateCorrectionPos();
  dev_t void calcTotMass();
	
  //__device__ vector_t & getVElem(const int &e, const int &n){return v[m_elnod[e*m_nodxelem+n]];}
  inline dev_t double  getVElem(const int &e, const int &n,const int &d){return v[3*m_elnod[m_nodxelem*e+n]+d];}  
  
  inline dev_t vector_t getV(const int &n){return make_vector_t(v[3*n], v[3*n+1], v[3*n+2]);}  
  
  
	void SolveChungHulbert();
  
  void SetDT(const double &dt_){dt=dt_;}
  void SetEndTime(const double &tf_){end_t=tf_;}
	int & getDim(){return m_dim;}
protected:
  double          m_tot_mass; //Only for testing, not needed
  dev_t int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
	int 						m_dim;
  int             m_node_count, m_elem_count;
	
	unsigned int 		*m_elnod;
  //unsigned int    *m_eloffset; //FROM FIRST ELEMENT NODE
	
	double 				 *x, *x_h; //Vector is double, x_h is used only on host
	double 				 *v; //CHANGED TO DOUBLE
	double 				 *a, *prev_a;
	double         *u, *u_dt;
  double         *v_h,*u_h, *a_h;
  double         *m_mglob, *m_mdiag, *m_ematm; //Global Matrix (TODO: MAKE SPARSE) and Diagonal masses, and element masses
  
	double 					*p, *rho, *rho_0;
  
  
  double          dt, end_t;
  
  double          m_alpha, m_beta, m_gamma;
  double          *vol_0, *vol;  //Element initial and current volume
  
  
  //TO NOT OCCUPY SO MUCH MEMORY
  int             *bcx_nod, *bcy_nod, *bcz_nod;
  double          *bcx_val, *bcy_val, *bcz_val;
  int             bc_count[3];
  
  ////// PARALLELIZATION
  int             Nproc;
  
  std::vector<double>     bcx_val_h,bcy_val_h,bcz_val_h;
  std::vector<int>        bcx_nod_h,bcy_nod_h,bcz_nod_h;
  
  Material_ **mat; //pointer to material of each particle
  Material_ *materials; //All materials 
  
  bool            m_red_int;  //Reduced integration, 1 GAUSS POINT
  int             m_gp_count; //Assuming constant gauss points
  int             m_nodxelem;
  
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
  
  
  ////////////////////////////// REDUCTION THINGS ////////////////////////////////////////
  //THESE ARE USED IF NODAL FORCES WANTED TO BE REDUCED (SUMMED) NOT BY ATOMIC (GPU) OR LOCK (OPENMP CPU)
  bool            m_nonlock_sum;
  int            *m_nodel;           // NODE ELEMENTS: ELEMENTS SHARED BY EACH NODE [nodxelem* node_count] call:  m_nodel[nod,elcount] =EL i,e m_nodel[nod_offset+elcount]
  int            *m_nodel_loc;        //
  int            *m_nodel_offset;    //OFFSET OF THE
  int            *m_nodel_count;    
	
	
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
__global__ void calcElemInitialVolKernel(Domain_d *dom_d);

__global__ void AssignMatAddressKernel(Domain_d *dom/*, Material_ *mat*/);

__global__ void calcTotMassKernel(Domain_d *dom_d);

__global__ void calcMassDiag(Domain_d *dom_d);

__global__ void calcAccelKernel(Domain_d *dom_d);

__global__ void printVecKernel(Domain_d *dom_d, double *);
__global__ void printSymmTensKernel(Domain_d *dom_d,  double *);

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
#endif
