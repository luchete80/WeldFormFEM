/*************************************************************************/
/*  ImpDomain_d_BORRAR.h                                         */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                              */
/*  https://www.opensourcemech.com                                                                */
/*                                                                       */
/*  Copyright (c) 2025-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/


#ifndef _IMPDOMAIN_D_CUH_
#define _IMPDOMAIN_D_CUH_

//////////////////// COMMON DOMAIN FOR CPU/GPU ///////////////////
#include "Domain_d.h"

#include <stdio.h>


#include "utils.h"

#include "Material.cuh"

#include <vector>
#include <string>
#include "lsdynaReader.h"


#define _QUA2D_ 0
#define _TRI2D_ 1
#define _TET2D_ 2
     

                                      
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


class Matrix;



namespace MetFEM{

class TriMesh_d;
class VTKWriter;
class VTUWriter;
class ReMesher;
class ReMesher;
class Solver;
class Solver_Eigen;

class ImpDomain_d :public Domain_d {
  friend class VTKWriter;
  friend class VTUWriter;
  friend class ReMesher;
  friend class Solver;
  friend class Solver_Eigen;
  
public:
  ImpDomain_d (std::string);
  ImpDomain_d (){
    m_axisymm_vol_weight = false;
    m_domtype = _Plane_Strain_;
    bc_count[0]=bc_count[1]=bc_count[2]=0;
    contact = false;
    m_thermal = false;
    m_remesh_interval = 1e10;
    
    //DEFAULTS
    m_dim = 3;
    m_gp_count = 1;
    m_nodxelem = 4;
  }
  //~ void setNproc(const int &n){Nproc=n;}
  virtual void SetDimension(const int &node_count, const int &elem_count); //ELEM TYPE???
  //~ void AddBoxLength(vector_t const & V, vector_t const & L, const double &r, const bool &red_int = true, const bool &tritetra = false);
  //~ void CreateFromLSDyna(LS_Dyna::lsdynaReader &reader);
  //~ ///// (CUDA HOST) FUNCTIONS 
  //~ inline double2 getPosVec2(const int &n){
    //~ return make_double2(x[m_dim*n], x[m_dim*n+1]);
    //~ };
  //~ inline vector_t getPosVec3  (const int &n){return make_vector_t(x[m_dim*n], x[m_dim*n+1], x[m_dim*n+2]);}; //the same
  

  //~ #ifdef CUDA_BUILD
  //~ inline vector_t getPosVec (int &n){return make_vector_t(x_h[m_dim*n], x_h[m_dim*n+1], x_h[m_dim*n+2]);};
  //~ inline double2  getPosVec2_h(const int &n){return make_double2 (x_h[m_dim*n], x_h[m_dim*n+1]);}; //the same
  //~ inline vector_t getPosVec3_h(const int &n){return make_vector_t(x_h[m_dim*n], x_h[m_dim*n+1], x_h[m_dim*n+2]);}; //the same
  //~ #else
  //~ inline vector_t getPosVec3_h(const int &n){return make_vector_t(x[m_dim*n], x[m_dim*n+1], x[m_dim*n+2]);}; //the same of getPosVec3
  //~ inline double2 getPosVec2_h(const int &n){
    //~ return make_double2(x[m_dim*n], x[m_dim*n+1]);
    //~ };
  //~ #endif

  //~ void setTriMesh(TriMesh_d *m){trimesh = m;}
  
  //~ void SearchExtNodes();

  //~ inline vector_t getAccVec (const int &n){return make_double3(a[m_dim*n], a[m_dim*n+1], a[m_dim*n+2]);};  
  //~ inline vector_t getContForceVec(const int &n){return make_vector_t(contforce[m_dim*n], contforce[m_dim*n+1], contforce[m_dim*n+2]);}
  
  //~ #ifdef CUDA_BUILD
  //~ inline vector_t getDispVec(const int &n){return make_vector_t(u_h[m_dim*n], u_h[m_dim*n+1], u_h[m_dim*n+2]);};
  //~ inline vector_t getVelVec (const int &n){return make_vector_t(v_h[m_dim*n], v_h[m_dim*n+1], v_h[m_dim*n+2]);};
  //~ inline vector_t getIntForceVec(const int &n){/*return make_vector_t(m_fi[m_dim*n], m_fi[m_dim*n+1], m_fi[m_dim*n+2]);*/};
  //~ #else
  //~ //inline vector_t getAccVec (const int &n){return make_vector_t(a[m_dim*n], a[m_dim*n+1], a[m_dim*n+2]);};
  //~ inline vector_t getVelVec (const int &n){return make_vector_t(v[m_dim*n], v[m_dim*n+1], v[m_dim*n+2]);};
  //~ inline vector_t getDispVec(const int &n){return make_vector_t(u[m_dim*n], u[m_dim*n+1], u[m_dim*n+2]);};
  //~ inline vector_t getIntForceVec(const int &n){return make_vector_t(m_fi[m_dim*n], m_fi[m_dim*n+1], m_fi[m_dim*n+2]);};
  //~ //inline vector_t getContForceVec(const int &n){return make_vector_t(contforce[m_dim*n], contforce[m_dim*n+1], contforce[m_dim*n+2]);}
  //~ #endif
  
  //~ void printVec(double*);
  //~ void printSymmTens( double *);
  
  //~ void setNodElem(int *elnod);//For assembly and parallel processing
  
  //~ ///// ALREADY ALLOCATED
  //~ void setNode(const int &i, const double &_x, const double &_y, const double &_z);
  
  //~ void WriteToVTK(char *);
  //~ int WriteToCSV(char *);
  
  //~ void calcElemJacobian ();
  //~ void calcElemJAndDerivatives/*_FullInt*/();
  //~ void calcElemMassMat();

	//~ void calcElemStrainRates();
  //~ void calcElemForces();
  //~ void calcElemHourglassForces();
  //~ void calcElemPressure(); //FROM STRAIN
  //~ void calcElemPressureFromJ();
  //~ void calcElemPressureANP(); //AVERAGE NODAL POINT

  //~ void assemblyMassMatrix();
  //~ void assemblyForces();
  
  //~ host_   void AddBCVelNode(const int &node, const int &dim, const double &val);
  //~ host_   void AllocateBCs();
  
  //~ void calcMassDiagFromElementNodes(const double &rho); // To use existing array, NOT IN PARALLEL

  //~ void ImposeBCA(const int dim); /// DO NOT USE REFERENCESSS!!!!!!
  //~ host_ void ImposeBCAAllDim();

  //~ void ImposeBCV(const int dim); /// DO NOT USE REFERENCESSS!!!!!!
  //~ host_ void ImposeBCVAllDim();
  
  //~ //void calcMinEdgeLength();
  
  //~ ///// ATENTION! THIS IS Deriv x DETJ
  //~ inline double & getDerivative(const int &e, const int &gp, const int &i, const int &j); //I AND J ARE: DIMENSION AND NODE
  //~ inline void     setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &); //I AND J ARE: DIMENSION AND NODE
  //~ inline void     setDerivative(const int &e, const int &gp, Matrix *); //I AND J ARE: DIMENSION AND NODE
  
  //~ void setDensity(const double &r);
  
  //~ // template <typename T>
  //~ // void setVector(,const double &r);
  
	//~ int threadsPerBlock, blocksPerGrid; //TO BE USED BY SOLVER
	
	//~ const int & getElemCount()const{return m_elem_count;}
  //~ unsigned int& getElemNode(const int &e, const int &n){return m_elnod[m_nodxelem*e+n];}
	//~ const int & getNodeCount()const{return m_node_count;}
  //~ vector_t getNodePos3(const int &n){
    //~ #ifdef CUDA_BUILD
    //~ return make_vector_t (x_h[m_dim*n],x_h[m_dim*n+1],x_h[m_dim*n+2]);
    //~ #else
    //~ return make_vector_t (x  [m_dim*n],x  [m_dim*n+1],x  [m_dim*n+2]);      
    //~ #endif
  //~ }


  //~ double2 getNodePos2(const int &n){
    //~ #ifdef CUDA_BUILD
    //~ return make_double2 (x_h[m_dim*n],x_h[m_dim*n+1]);
    //~ #else
    //~ return make_double2 (x  [m_dim*n],x  [m_dim*n+1]);      
    //~ #endif
  //~ }

  
  //~ void AssignMaterial (Material_ *material_h); //Create and copy material
  
  void AssignMatAddress(); //Assign particle data to material array to zero array
  void AssignMatAddress_(); //Assign particle data to material array to zero array
  
  //~ void CalcStressStrain(double dt);
  void Calc_Elastic_Stress (const double dt);
  
  //~ ///// AXISYMM VARS
  //~ void Calc_Element_Radius(); //For axisymm
  //~ dom_type m_domtype;
  //~ bool m_axisymm_vol_weight;
  //~ void setDomtype();

  //~ //TODO: CHANGE THIS
  //~ inline double & getSigma  (const int e, const int gp, int i, int j){int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}};/*printf("i j symm idx %d %d %d\n",i,j,symm_idx[i][j]);*/if (j<3) return m_sigma   [e*m_gp_count*6 + symm_idx[i][j]];}
  //~ inline double & getStrRate(const int e, const int gp, int i, int j){int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}}; if (j<3) return m_str_rate[e*m_gp_count*6 + symm_idx[i][j]];}
  
  //~ void CalcElemInitialVol();
  
  //~ void InitUVA(); 
  
  //~ void CalcElemVol();
  //~ //Specially for tetra
  //~ void CalcNodalVol();
  //~ void CalcNodalMassFromVol();
  
  //~ void calcElemDensity();
    
  //~ void calcAccel();
  //~ void UpdatePrediction();
  //~ void UpdateCorrectionAccVel();
  //~ void UpdateCorrectionPos();
  //~ void calcTotMass();
  //~ void InitElemValues(double *arr, double val = 0.0);
  //~ void InitStressesFromMat();
    
  //~ inline double getRadius(int e, int gp){return m_radius[e*m_gp_count+gp];}
	
  //~ //__device__ vector_t & getVElem(const int &e, const int &n){return v[m_elnod[e*m_nodxelem+n]];}
  //~ inline double  getVElem(const int &e, const int &n,const int &d){return v[m_dim*m_elnod[m_nodxelem*e+n]+d];}  
  
  //~ inline vector_t getV(const int &n){return make_vector_t(v[m_dim*n], v[m_dim*n+1], v[m_dim*n+2]);}  
  
  
	void Solve();
	 void ElasticSolve();
  //~ void setdtOut(const double &t){m_dtout=t;}
  //~ void SetDT(const double &dt_){dt=dt_;}
  //~ void SetEndTime(const double &tf_){end_t=tf_;}
	//~ int & getDim(){return m_dim;}
  //~ void setProcCount(const int &np){Nproc = np;} //for CPU only
  
  //~ void CalcContactForces(); //standard, not predicted
  //~ void setContactOn(){contact = true;}
  //~ bool isContactOn(){return contact;}
  //~ TriMesh_d* getTriMesh(){return trimesh;}
  //~ void InitValues();
  //~ double getMinLength(){return m_min_length;}
  

  //~ bool m_auto_contact; // if not define external nodes
  //~ void setRemeshInterval(int i) {m_remesh_interval = i;}
  //~ void Free();
  // IMPLICIT FUNCTIONS
  //--------------------------------------------------------------------------------------------------------------------------------
  Matrix getElemBMatrix(const int &e);
  void CalcMaterialStiffElementMatrix();
  void CalcGeomStiffElementMatrix();  
  void CalcElemIntForces();

  //--------------------------------------------------------------------------------------------------------------------------------

  
protected:
  double          m_tot_mass; //Only for testing, not needed
  int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
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
  
  //////Thermal
  bool m_thermal;
  double *T;   /// temperature
  double *m_dTedt; //elem node
  double m_min_length;
  double m_min_height;
  double m_cfl_factor;
  int             m_remesh_interval;

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
  
  //////////////////////////////// IMPLICIT THINGS
  Matrix **m_Kmat;   //MATERIAL PART OF 
  Matrix **m_Kgeo;   //MATERIAL PART OF   
  
  //// DEFORMATION GRADIENT
  Matrix **m_Fel;    //def
  
  
  Solver* m_solver;
  
};




}; //Namespace

//#include "Mesh.h"
//#include "Contact.C"

#endif
