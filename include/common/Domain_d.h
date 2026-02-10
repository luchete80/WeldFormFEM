/*************************************************************************/
/*  Domain_d.h                                                   */
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

#include "Tensor.h"

 
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

#define avgScalarDom(v,a,dim, dom)    for (int n=0;n<dom->m_node_count;n++){\
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

enum class TimeInt {
    IMPLICIT,       // Implicit time integration
    EXPLICIT,       // Explicit time integration
    QUASI_STATIC,   // No inertial effects
    DYNAMIC         // Generic dynamic (could be either)
};

struct StabilizationParams {
    double alpha_free;
    double alpha_contact;
    double hg_coeff_free;
    double hg_coeff_contact;
    double av_coeff_div;
    double av_coeff_bulk;
    double log_factor;
    double pspg_scale;
    double p_pspg_bulkfac;
    double J_min;
    double hg_visc;
    double hg_stiff;
};

enum class PlasticityType {
    Hardening,     // Elasto-plástico con endurecimiento
    Perzyna,       // Rígido-plástico viscoplástico
    Norton         // Similar a Perzyna, otro flujo viscoplástico
};
 

// Structure to define a face with 4 nodes
//FOR HEXA IS 4
//PARALLELIZE WITH GPU : TODO
//// ORIGINALLY MEANS HEXA

// //OLD FOT HEXA, CHANGE IT
#define MAX_FACE_NODES 4


// Triángulo 2D: 3 aristas de 2 nodos
/*__constant__*/
static const int tri_edges[3][MAX_FACE_NODES] = {
    {0,1,-1,-1}, {1,2,-1,-1}, {2,0,-1,-1}
};

// Quad 2D: 4 aristas de 2 nodos
/*__constant__*/
static int quad_edges[4][MAX_FACE_NODES] = {
    {0,1,-1,-1}, {1,2,-1,-1}, {2,3,-1,-1}, {3,0,-1,-1}
};

// Tetra 3D: 4 caras de 3 nodos
/*__constant__*/
static int tetra_faces[4][MAX_FACE_NODES] = {
    {0,1,2,-1}, {0,1,3,-1}, {1,2,3,-1}, {0,2,3,-1}
};

// Hexa 3D: 6 caras de 4 nodos
/*__constant__*/
static int hexa_faces[6][MAX_FACE_NODES] = {
    {0,1,2,3}, {4,5,6,7}, {0,1,5,4},
    {2,3,7,6}, {0,3,7,4}, {1,2,6,5}
};

struct Face {
    int nodes[MAX_FACE_NODES];  // Fast to avoid dynamic assignment
    int n_nodes;                // How manynodes are used
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
  double3 start, end;
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
    Nproc = 1;    //USED ON IMPLICIT
    m_axisymm_vol_weight = false; // DO NOT CHANGE, SPECIALLY FOR IMPLICIT SOLVER.
    m_domtype = _3D_;
    m_dim = 3;
    //IF CONTACT 
    ELNOD=4;   //ORIGINALLY 8
    FACENOD=3;  //ORIGINALLY 4
    ELFAC=4;   //ORIGINALLY 6
    local_faces = tetra_faces; 


    bc_count[0]=bc_count[1]=bc_count[2]=0;
    contact = false;
    m_thermal = false;
    m_remesh_interval = 1e10;
    m_contPF = 0.1;
    m_fixed_dt = false;
    m_press_algorithm = 0;
    
    m_remesh_min_fr = 0.4;
    m_remesh_max_fr = 2.0;    

    m_remesh_min_pl_strain = 1.0e-1;
    m_remesh_max_pl_strain = 1.0e6;
    m_remesh_map_vel = false;
    m_remesh_map_acc = false;
    m_remesh_max_count  =1e6;
    
    m_remesh_damp_vel = 0.02;
    
    m_plheatfraction = 0.9;
  
    //~ m_stab.alpha_free       = 0.3;
    //~ m_stab.alpha_contact    = 0.7;
    //~ m_stab.hg_coeff_free    = 0.2;
    //~ m_stab.hg_coeff_contact = 0.08;
    //~ m_stab.av_coeff_div     = 0.15;
    //~ m_stab.av_coeff_bulk    = 0.15;
    //~ m_stab.log_factor       = 0.8;
    //~ m_stab.pspg_scale       = 0.2;
    //~ m_stab.p_pspg_bulkfac   = 0.05;

    m_stab.alpha_free       = 0.0;
    m_stab.alpha_contact    = 0.0;
    m_stab.hg_coeff_free    = 0.0;
    m_stab.hg_coeff_contact = 0.0;
    m_stab.av_coeff_div     = 0.0;
    m_stab.av_coeff_bulk    = 0.0;
    m_stab.log_factor       = 0.0;
    m_stab.pspg_scale       = 0.0;
    m_stab.p_pspg_bulkfac   = 0.0;
    m_stab.J_min            = 0.0;
    m_stab.hg_visc        = 0.0;
    m_stab.hg_stiff       = 0.1;
    
    if (m_dim == 2 && m_nodxelem == 4) m_stab.hg_visc = 0.1;
    
    abs_bc_initialized = false;
    
    m_symm[0]=m_symm[1]=m_symm[2]=false;
        
  }

    // ==================== REMESH ====================
    struct FilterFlags {
        bool enable_warmup_damping;      // Damping durante warmup
        bool enable_velocity_peaks_correction; // Corrección de picos de velocidad
        bool enable_stress_filtering;    // Filtrado de tensiones
        bool enable_energy_based_limiter; // Limitador basado en energía
        bool enable_post_remesh_filter;  // Filtro post-remallado
        bool enable_strain_rate_smoothing; // Suavizado de tasas de deformación
        bool enable_laplacian_smoothing; // Suavizado laplaciano de campos
        bool enable_adaptive_dt_recovery; // Recuperación adaptativa de dt
        
        // Valores por defecto (conservadores)
        FilterFlags() : 
            enable_warmup_damping(true),
            enable_velocity_peaks_correction(true),
            enable_stress_filtering(true),
            enable_energy_based_limiter(true),
            enable_post_remesh_filter(true),
            enable_strain_rate_smoothing(false), // Más agresivo
            enable_laplacian_smoothing(false),   // Más agresivo  
            enable_adaptive_dt_recovery(true)
        {}
    };
    
    struct FilterParameters {
        // Warmup parameters
        int warmup_steps;                // Pasos de warmup después de remallado
        double warmup_damping_factor;    // Factor de damping durante warmup (0-1)
        
        // Velocity control
        double max_allowed_velocity;     // Velocidad máxima permitida
        double velocity_correction_aggressiveness; // 0-1: agresividad de corrección
        
        // Stress filtering  
        double stress_smoothing_factor;  // Factor de suavizado de tensiones (0-1)
        int stress_filter_interval;      // Cada cuántos pasos filtrar tensiones
        
        // Energy-based limits
        double max_energy_increase_ratio; // Máximo aumento permitido de energía
        double energy_correction_strength; // Fuerza de corrección energética
        
        // Post-remesh specific
        double post_remesh_blend_duration; // Duración blending post-remallado (0-1)
        double plastic_strain_threshold;  // Umbral de deformación plástica para filtrado
        
        int trans_step_count;
        
        // Default values
        FilterParameters() :
            warmup_steps(10),
            warmup_damping_factor(0.8),
            max_allowed_velocity(50.0),
            velocity_correction_aggressiveness(0.7),
            stress_smoothing_factor(0.3),
            stress_filter_interval(5),
            max_energy_increase_ratio(2.0),
            energy_correction_strength(0.5),
            post_remesh_blend_duration(0.1),
            plastic_strain_threshold(0.05),
            trans_step_count(50)
        {}
    };
    
    // Miembros en Domain_d
    FilterFlags m_filter_flags;
    FilterParameters m_filter_params;
    
    // ==================== MÉTODOS DE FILTRADO MODULARES ====================
    void applyWarmupDamping(double s_wup);
    void correctVelocityPeaks();
    void filterStresses();
    void applyEnergyBasedLimiting(double dEkin, double dEint);
    void applyPostRemeshFiltering();
    void smoothStrainRates();
    void applyLaplacianSmoothing();
    void adaptiveDTRecovery();
    
  ////////////// REMESH

  TimeInt m_timeint_type;
  void setNproc(const int &n){Nproc=n;}
  void SetDimension(const int &node_count, const int &elem_count); //Common
  void SetDimensionExplicit(const int &node_count, const int &elem_count); //ELEM TYPE???
  void SetDimensionImplicit(const int &node_count, const int &elem_count); //TO NOT USE virtual (GPU issues)
  void AddBoxLength(vector_t const & V, vector_t const & L, const double &r, const bool &red_int = true, const bool &tritetra = false);
  void CreateFromLSDyna(LS_Dyna::lsdynaReader &reader);
  void setSymm(const int &d){m_symm[d]=true;}
  ///// (CUDA HOST) FUNCTIONS 
  inline dev_t double2 getPosVec2(const int &n){
    return make_double2(x[m_dim*n], x[m_dim*n+1]);
    };
  inline dev_t vector_t getPosVec3  (const int &n){
    double z = 0.0;
    if (m_dim == 3) z = x[m_dim*n+2];
    return make_vector_t(x[m_dim*n], x[m_dim*n+1], z);}; //the same
  

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

  inline vector_t getAccVec (const int &n){ return make_double3(a[m_dim*n], a[m_dim*n+1], m_dim==3?a[m_dim*n+2]:0.0);}  
  inline vector_t getContForceVec(const int &n){return make_vector_t(contforce[m_dim*n], contforce[m_dim*n+1], contforce[m_dim*n+2]);}

  void setContForceVec(const int &n, const int &d, const double &v){contforce[m_dim*n+d]=v;}
  
  #ifdef CUDA_BUILD
  inline vector_t getDispVec(const int &n){return make_vector_t(u_h[m_dim*n], u_h[m_dim*n+1], u_h[m_dim*n+2]);};
  inline vector_t getVelVec (const int &n){return make_vector_t(v_h[m_dim*n], v_h[m_dim*n+1], v_h[m_dim*n+2]);};
  inline vector_t getIntForceVec(const int &n){/*return make_vector_t(m_fi[m_dim*n], m_fi[m_dim*n+1], m_fi[m_dim*n+2]);*/};
  #else
  //inline vector_t getAccVec (const int &n){return make_vector_t(a[m_dim*n], a[m_dim*n+1], a[m_dim*n+2]);};
  inline vector_t getVelVec (const int &n){return make_vector_t(v[m_dim*n], v[m_dim*n+1], m_dim==3?v[m_dim*n+2]:0.0);};
  inline vector_t getDispVec(const int &n){
    double z = 0.0;
    if (m_dim == 3)
      z = u[m_dim*n+2];
    return make_vector_t(u[m_dim*n], u[m_dim*n+1], z);
  }
  inline vector_t getIntForceVec(const int &n){return make_vector_t(m_fi[m_dim*n], m_fi[m_dim*n+1], m_dim==3?m_fi[m_dim*n+2]:0.0);};
  //inline vector_t getContForceVec(const int &n){return make_vector_t(contforce[m_dim*n], contforce[m_dim*n+1], contforce[m_dim*n+2]);}
  #endif
  
  dev_t void printVec(double*);
  dev_t void printSymmTens( double *);
  
  void setNodElem(int *elnod);//For assembly and parallel processing
  
  void WriteToVTK(char *);
  int WriteToCSV(char *);
  
  dev_t void calcElemJacobian ();
  dev_t void calcElemJAndDerivatives/*_FullInt*/();
  dev_t void calcElemJAndDerivatives_Tet_SSRI();
  dev_t void calcElemMassMat();

	dev_t void calcElemStrainRates();
  dev_t void calcElemForces();
  dev_t void calcElemHourglassForces();
  dev_t void calcElemPressure_Hex(); //FROM STRAIN
  dev_t void calcElemPressure(); //FROM STRAIN
  dev_t void calcElemPressureRigid(const double &K);
  dev_t void calcElemElasticJ();
  dev_t void calcElemPressureLocal(); //SF Style
  dev_t void calcElemPressureCaylak();
  dev_t void smoothPressureLaplacian();
  
  dev_t void smoothFieldLaplacian(double *, int dim = 3);
  
  int AddBCVelZone(const vector_t &start, const vector_t &end, const vector_t &vel);
    
  dev_t void smoothPressureField(double gamma); //Laplace Smooth
  dev_t void calcElemPressureFromJ();
  dev_t void calcElemPressureANP(); //AVERAGE NODAL POINT
  dev_t void calcElemPressureElementBased();
  dev_t void calcElemPressureANP_Nodal();
  dev_t void calcElemPressureANP_Nodal_HG();
  dev_t void calcElemPressureANP_Nodal_Stab();
  
  dev_t void calcElemStrGradF();
  
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
    return make_vector_t (x  [m_dim*n],x  [m_dim*n+1],m_dim==3?x  [m_dim*n+2]:0.0);      
    #endif
  
  }


  double2 getNodePos2(const int &n){
    #ifdef CUDA_BUILD
    return make_double2 (x_h[m_dim*n],x_h[m_dim*n+1]);
    #else
    return make_double2 (x  [m_dim*n],x  [m_dim*n+1]);      
    #endif
  }

  void setFixSymm(); //Set Nodes for fixing 
  void AssignMaterial (Material_ *material_h); //Create and copy material
  
  dev_t void AssignMatAddress(); //Assign particle data to material array to zero array
  void AssignMatAddress_(); //Assign particle data to material array to zero array
  
  dev_t void smoothDevStrainRates(double beta);
  dev_t void CalcStressStrain(double dt);
  dev_t void CalcStressStrainRigidViscoPlastic(double dt);
  dev_t void CalcStressStrainElastoViscoPlastic(double dt);
  
  dev_t Matrix CalcElementStressAndTangent(int e, double dt); ///// PER ELEMENT (FOR IMPLICIT)
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
  
  inline dev_t vector_t getV(const int &n){return make_vector_t(v[m_dim*n], v[m_dim*n+1], m_dim==3?v[m_dim*n+2]:0.0);}  
  
  dev_t void calcThermalExpansion();
  
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

  inline double axisymWeight(double r_gp) {
      return (m_domtype == _Axi_Symm_) ? 2.0 * M_PI * r_gp : 1.0;
  }
  inline double computeDV(double detJ, double wgp, double r_gp) {
      double dV = detJ * wgp;
      if (m_domtype == _Axi_Symm_) {
          dV *= 2.0 * M_PI * r_gp;
      }
      return dV;
  }
  void setAxiSymm(const bool &vol_weight = false){
    m_dim = 2;
    m_domtype = _Axi_Symm_;
    m_axisymm_vol_weight = vol_weight;
  }
  StabilizationParams m_stab;
  int m_press_algorithm;
  
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
  
  dev_t void CalcExtFaceAreas();
  std::ofstream out_file;
  void setCFL(const double &f){m_cfl_factor = f;}
  
    // IMPLICIT FUNCTIONS
  //--------------------------------------------------------------------------------------------------------------------------------
  Matrix getElemBMatrix(const int &e);
  void CalcMaterialStiffElementMatrix();
  void CalcGeomStiffElementMatrix();  
  void CalcElemIntForces();
	void ElasticIncSolve(); ///TESTING INCREMENTAL ELASTIC ALL IN ONE FUNCTION
	void ElasticSolve();
  
  void SolveImplicitGlobalMatrix();
  /// DISPLACEMENT (STATIC) FUNCTIONS
  void SolveStaticQS_EPHV(); ///// ELASTOPLASTIC HARDENING, VIA VELOCITIES
  void SolveStaticQS_UP();
  void SolveStaticQS_V();
  
  void Solve_Martins_Picard();
  
  
  void SolveElastic(); //SIMPLE 
  
  ///// SIMILAR, INCREMENTAL BCs which Correct prescribed BCs for incremental solver.
  void CalcIncBCU(int dim);
  void CalcIncBCV(int dim/*, double load_factor*/);
  Matrix getConsistentPlasticTangentMatrix(const tensor3& s_trial, double sig_trial, 
                                                  double G, double H);

  dev_t void calcElemForces(const int &e); //FOR IMPLICIT
  
  inline void getElementVelocityDOFs(int e, std::vector<int>& dofs) {
    const int ndof = m_nodxelem * m_dim;
    dofs.resize(ndof);

    for (int a = 0; a < m_nodxelem; ++a) {
        const int node = getElemNode(e, a);
        for (int d = 0; d < m_dim; ++d) {
            dofs[a*m_dim + d] = node*m_dim + d;
        }
    }
  }
  inline int getElementPressureDOF(int e) {
    return m_dim*m_node_count + e;
  }

  //--------------------------------------------------------------------------------------------------------------------------------
  //IF USING INCREMENTAL BCs
  void setOriginalBCs() {
                printf("Ckecking\n");
      if (!abs_bc_initialized) {
          printf("Inititalizing original bcs\n");
          // Guardar copia de los valores originales
          original_bcx_val.assign(bcx_val, bcx_val + bc_count[0]);
          original_bcy_val.assign(bcy_val, bcy_val + bc_count[1]);
          original_bcz_val.assign(bcz_val, bcz_val + bc_count[2]);
          abs_bc_initialized = true;
      }
  }
  
  void setU(const int &n, const int &d, const double &val){u[m_dim*n+d] = val;}

  
  void setFixedDt(const bool &f){m_fixed_dt = f;}
  const bool &ifFixedDt()const{return m_fixed_dt;}
  double m_artifvisc[2];
  void setContactPF(const double &pf){m_contPF = pf;}

  __device__ void ApplyGlobalSprings();
  void ApplyGlobalDamping(double damping_factor,double vel_ref = 0.5);
  
  double  m_remesh_min_pl_strain,m_remesh_max_pl_strain;
  bool m_remesh_map_vel,m_remesh_map_acc;
  int m_remesh_max_count;
  double m_remesh_damp_vel;
  
  double m_remesh_min_fr,m_remesh_max_fr; 
  
  double m_plheatfraction;
  
  void set2DFacesValues(){      
    ELNOD   = 4;   //SAME AS 3D
    FACENOD = 2;   //4 in 3D tetra
    ELFAC   = 4;   //Quad faces
      //dom_d->local_faces = tetra_faces; 
    local_faces = quad_edges;
  }
  
  ///// REMESH
  dev_t void BlendStresses(const double &s, const double &pl_strain_max);
  dev_t void BlendField(const double &s, const int size, const int &d, double *prev, double *curr);
  dev_t void postRemeshGlobFilter();
  dev_t void SmoothDeviatoricStress(double alpha);
  dev_t void computeEnergies(double dt,
                       //const double* f_visc,  // optional nodal viscous forces [3*node_count] or nullptr
                       double &Ekin, double &Eint, double &Evisc);
                       
  dev_t double getPtrMax(double *v, const int &size, const int &dim);
  
  double3 ComputeCentroid() const; 
  
  int getBCCount(const int &d){if (d==0)return bcx_nod_h.size(); else if(d==1)return bcy_nod_h.size(); else if (d==2)return bcz_nod_h.size(); else return 0;}
  
  void CorrectLocalVelocityPeaks();
  
  void setName(std::string name){m_name = name;}
  std::string getName(){return m_name;}
  
  
  
  double m_remesh_min_frac  = 0.4;    // min_size = 0.35*h0
  double m_remesh_max_frac  = 2.0;     // max_size = 2.0*h0
  double m_remesh_eps_ref   = 1.0;    // plastic strain que consideras "alto" (ajustar)
  double m_remesh_beta      = 4.0;     // controla rapidez del refinamiento exponencial
  int    m_remesh_type      = 0;        // LINEAR (HYPERBOLIC)
  
  bool m_devElastic = true;
  PlasticityType m_plastType = PlasticityType::Hardening; //PlasticityType::Perzyna //Norton
  
  double m_max_edot = 1.0e6;
  double m_min_edot = 1.0e-10;
  double m_eta_visc = 1.0e-3;
  double m_m_visc   = 1.0;
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
protected:
  std::string     m_name;
  double          m_tot_mass; //Only for testing, not needed
  dev_t int symm_idx[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
	int 						m_dim;
  int             m_node_count, m_elem_count;
  bool            m_symm[3];
	
	unsigned int 		*m_elnod;              /// TO SHIFT 
  
  
  unsigned int   *m_elnodoffset; //FROM FIRST ELEMENT NODE
	
	double 				 *x, *x_h; //Vector is double, x_h is used only on host
	double 				 *v, *prev_v; //PREV_V IS ONLY USED ON IMPLICIT ANALYSIS
	double 				 *a, *prev_a;
	double         *u, *u_dt;
  double         *v_h,*u_h, *a_h;
  double         *m_mglob, *m_mdiag, *m_ematm; //Global Matrix (TODO: MAKE SPARSE) and Diagonal masses, and element masses
  
	double 					*p, *rho, *rho_0;
  
  double          *ut_prev;   //FOR HISTORICAL CONTACT
  
  double          *m_vprev;
  
  //////Thermal
  bool m_thermal;
  double *T;   /// temperature
  double *m_dTedt; //elem node
  double *ps_energy;
  double *q_cont_conv;
  double *node_area; //THERMAL CONTACT
  double  contact_hc; //MOVETO CONTACT INTERFACE
  double *m_elem_length; //FOR CONTACT
  
  double *m_elem_min_angle,*m_elem_max_angle;
  double m_min_Jnorm;
  int  m_bad_elem_count;
  
  double *p_node; //FOR CONTACT AND ANP (ANP will be overriden by HG later)
  double *m_elem_area;
  
  int 	 *m_mesh_in_contact; //ID (Pos) with mesh in contact
  
  double m_min_length;
  double m_min_height;
  
  double m_min_angle,m_max_angle;
  double m_cfl_factor;
  int             m_remesh_interval;
  double          m_remesh_length;

    

  double          *pl_strain, *sigma_y;
  
  
  double          dt, end_t;
  
  double          m_alpha, m_beta, m_gamma;
  double          *vol_0, *vol;  //Element initial and current volume
  
  double          *m_voln,*m_vol_0n; //Nodal values, useful for TETRA nodal Average Nodal Pressure
  double          m_vol_tot;
  
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
  double        m_pl_energy;
  
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
  double          *m_sigma, *m_tau, *m_eps;
	double          *m_radius;
  
  
  ///// FOR PLASTIC/ELASTIC DECOMPOSITION
  double          *m_Jel; // Almacenar J elástico por elemento
  double          *m_Jpl;  // J plástico por elemento (para cálculo nodal)
  double          *m_F;
  double          *m_Fp;
  
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
  
  int ELNOD;   // nodos por elemento
  int FACENOD; // nodos por cara
  int ELFAC;   // caras por elemento
  const int (*local_faces)[MAX_FACE_NODES];
    
  /////// IMPLICIT
  //////////////////////////////// IMPLICIT THINGS
  Matrix **m_Kmat;   //MATERIAL PART OF 
  Matrix **m_Kgeo;   //MATERIAL PART OF   
  
  //// DEFORMATION GRADIENT
  Matrix **m_Fel;    //def
  
  
  Solver* m_solver;
  double *x_old;
  double m_contPF = 0.2;


  /////////////////////////////// ROLLBACK
      // Variables para guardar el estado anterior
  double          *m_sigma_prev, *pl_strain_prev,*m_tau_prev;
  double          *m_str_rate_prev;
  
  double* m_u_prev;         // Desplazamientos
  double* m_v_prev;         // Velocidades
  double* m_a_prev;         // Aceleraciones
  double* m_prev_a_prev;    // Aceleraciones del paso previo

  double* m_pl_strain_prev; // Deformación plástica
  double* m_vol_prev;       // Volúmenes elementales
  double* m_rho_prev;       // Densidades
  double* m_p_prev;         // Presiones

  //////////////////////////////////////////////////////////

  // Only used for incremental
  std::vector<double> original_bcx_val;
  std::vector<double> original_bcy_val; 
  std::vector<double> original_bcz_val;
  bool abs_bc_initialized;
  
  bool m_sym_axis[3]; //IN CASE OF PLANE SYMM, FOR REMESH.
  
  double m_dt_gap_min; //CONTACT
  double m_Kpen;
  bool m_autoKpen;
  
  double*  m_hg_q;



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
