/*************************************************************************/
/*  Solver_GlobalMatrix.C                                        */
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


#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

#include "Mesh.h"
#include "WallTimer.h"
#include "Matrix.h"

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>  // ofstream
#include <iomanip>

#ifdef BUILD_REMESH
#include "ReMesher.h"
#endif

#include <omp.h>

#include "Solver_Eigen.h"
#include "Tensor3.C"

using namespace std;


namespace MetFEM{






void Domain_d::CalcIncBCU(int dim/*, double load_factor*/) {
    int* bc_nodes = nullptr;
    double* current_bc_values = nullptr;
    std::vector<double>* original_values = nullptr;
    
    // Seleccionar arrays según dimensión
    switch (dim) {
        case 0: 
            bc_nodes = bcx_nod; 
            current_bc_values = bcx_val;
            original_values = &original_bcx_val;
            break;
        case 1: 
            bc_nodes = bcy_nod; 
            current_bc_values = bcy_val;
            original_values = &original_bcy_val;
            break;
        case 2: 
            bc_nodes = bcz_nod; 
            current_bc_values = bcz_val;
            original_values = &original_bcz_val;
            break;
        default: return;
    }
    
    // Aplicar BCs de forma incremental según load_factor
    par_loop (i, bc_count[dim]) {
        double target_disp = (*original_values)[i] /** load_factor*/;
        double current_disp = u[bc_nodes[i] * m_dim + dim];
        double delta_disp = target_disp - current_disp;
        
        //cout << "Writing bc "<< i <<", node "<<bc_nodes[i]<<endl;
        // Sobrescribir el valor en el array de BCs
        current_bc_values[i] = delta_disp; // ← Ahora guarda el incremento

    }
}


void host_ Domain_d::SolveStaticQS_UP(){
  WallTimer timer;

  std::ofstream of("Contact_Forces.csv", std::ios::out);
  
  int N;
	N = getElemCount();
  #if CUDA_BUILD
	threadsPerBlock = 256; //Or BlockSize
	//threadsPerBlock = 1; //Or BlockSize
	blocksPerGrid =				// Or gridsize
	(N + threadsPerBlock - 1) / threadsPerBlock;
	cout << "Blocks per grid"<<blocksPerGrid<<", Threads per block"<< threadsPerBlock<<endl;

  ////// MATERIAL
  cout << "Assignin material.."<<endl;
  AssignMatAddressKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
  cudaDeviceSynchronize();
  #else 
  AssignMatAddress();
  #endif


  
  cout << "Initializing Values ..."<<endl;
  cout << "USING STATIC SOLVER"<<endl;
  InitValues();
  cout << "Done. "<<endl;  
  double3 *m_v_orig;
  if (contact){
    ////TEMP; SMOOTH LOAD(ASSUMING CONTSTANT)
    m_v_orig = new double3 [trimesh->nodecount];
    for (int n=0;n<trimesh->nodecount;n++){
      m_v_orig[n]= trimesh->node_v[n];
    }
  }
  
  //cout << "Imposing BCVs.. "<<endl;
  for (int d=0;d<m_dim;d++){
    
    #ifdef CUDA_BUILD
    ////REMAINS TO INIT VELOCITIES
    N = bc_count[d];
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    ImposeBCVKernel<<<blocksPerGrid,threadsPerBlock >>>(this, d);
    cudaDeviceSynchronize();
    #else
      for (int n=0;n<m_node_count*m_dim;n++){
        v[n]=a[n]=u[n]=0.0;
      }
       ImposeBCV(d);
    #endif
  }

  cout << "Done."<<endl;
  
  ostringstream oss_out;
  
  

  calcElemJAndDerivatives();

  CalcElemInitialVol(); //ALSO CALC VOL
  
  CalcElemVol();
  printf("calc dens\n");
  calcElemDensity();
  cout <<"Done."<<endl;
  // if (m_dim == 3 && m_nodxelem ==4){
  // //Replaces PREVIOUS, INSTEAD MASS APPROACH, BUT STILL TO WORK FOR HEXAS
  // cout << "Calc tetra vol"<<endl;
  cout <<"Calc Nodal Volume"<<endl;
  CalcNodalVol(); //To calc nodal mass
  CalcNodalMassFromVol(); //Repla
  
  cout << "Done."<<endl;
  

  Time = 0.0;
  int step_count = 0;
  double tout = 0;
  
  bool remesh_ = false;
  int remesh_count = 0;
  const double RAMP_FRACTION = 1.0e-2;  // 0.1% of total time instead of 1%
  of << "t,f,area"<<endl;
  

  setOriginalBCs(); //For incremental
  
  ///// SOLVER THINGS.
  Solver_Eigen *solver = new Solver_Eigen();
  m_solver = solver;
  m_solver->setDomain(this);
  m_solver->Allocate();

  double dt_initial = end_t / 1000.0; // Initial guess
  double dt_min = end_t / 10000.0;   // Minimum allowable
  double dt_max = end_t / 2.0;      // Maximum allowable
  double dt = dt_initial;
  
  
  double *delta_v, *x_initial;
  
  #ifndef BUILD_GPU
    delta_v   = new double [m_dim * m_node_count];
    x_initial = new double [m_dim * m_node_count];
  #else
  
  #endif
  
  double tau_old[m_elem_count*6];
  double sig_old[m_elem_count*6];
  double pls_old[m_elem_count  ];


  double prev_x[m_node_count * m_dim];
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "------------------------MAIN LOOP---------------------------"<<endl;
  cout << "Time: "<<Time <<", End Time "<< end_t<<endl;
  cout << "Time Step: "<<dt <<endl;
  
  bool end_all = false;
  while (Time < end_t && !end_all) {

  for (int i=0;i<m_elem_count*6;i++){
    tau_old[i] = m_tau[i];
    sig_old[i] = m_sigma[i];
  }
  for (int i=0;i<m_elem_count;i++) pls_old[i] = pl_strain[i];
  
  ////// OR TIME
  if (step_count % 100 == 0){
    printf("Step %d, Time %f, End Time: %.4e, Step Time %.4e\n",step_count, Time, end_t, dt);  
    timer.click();
    //std::cout << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    std::cout << "Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    //std::cout << "CPU Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  }

  
  //cout << "Storing previous values"<<endl;
  memcpy(prev_v,    v, sizeof(double) * m_node_count * m_dim);

  memcpy(x_initial, x, sizeof(double) * m_node_count * m_dim);

  memcpy(prev_x, x, sizeof(double) * m_node_count * m_dim);
  //cout << "Done."<<endl;
  
  //cout << "Calc External Faces"<<endl;
  if (step_count % 10 == 0){
    //cout << "Calc ExtFace Areas"<<endl;
    CalcExtFaceAreas();
    //cout << "Done"<<endl;
    }



  // Newton-Raphson loop
  double tolerance = 1e-4; //dv tol
  double ftol = 1e-4;

  //~ double tol_force = 1e-3;    // 0.1% error en fuerzas  
  //~ double tol_disp = 1e-4;     // Desplazamientos
  //~ double tol_energy = 1e-4;   // Energía residual


  int max_iter = 20;

  double force_factor = 1.0e-3;//TO AVOID ILL CONDITIONING
  
  
  double prev_Rnorm;
  double alpha_damp= 1.0;

  ////delta_v: Pure NR correction term
  for (int i=0;i<m_node_count*m_dim;i++)delta_v[i]=0.0;
 
  
  double flat_fold[6*m_elem_count];
  
  int nr_iterations = 0;
  int conv_iter = 0;


  // Al inicio del paso de tiempo, guardar TODO el estado
  double x_old[m_node_count * m_dim];
  double u_old[m_node_count * m_dim]; 
  double v_old[m_node_count * m_dim];


  // Guardar al inicio del paso de tiempo
  memcpy(x_old, x, sizeof(double) * m_node_count * m_dim);
  memcpy(u_old, u, sizeof(double) * m_node_count * m_dim);
  memcpy(v_old, v, sizeof(double) * m_node_count * m_dim);

  
  double u_inc[m_node_count * m_dim]; 
  
  ////////////////////////////////////////////////////////////
  ////////////////////////// NR LOOP /////////////////////////
  bool converged = false;
  int iter = 0;
  bool end = false;
  while (!converged && !end){

    nr_iterations++;
    cout <<"ITER "<<iter<<endl;
    printf("Step %d, Time %f, End Time: %.4e, Step Time %.4e\n",step_count, Time, end_t, dt);  

    const double gamma = 0.5;   // Newmark parameter
    // (1) Update velocities (v = prev_v + delta_v)

  double maxv[]={0.0,0.0,0.0};
 
    for (int i = 0; i < m_node_count * m_dim; i++) {
        v[i] = prev_v[i] + delta_v[i]; ///v: Current total velocity 
    }
 
   for (int e=0;e<m_node_count;e++)
      for (int d=0;d<3;d++)
            if (abs(v[e*m_dim+d])>maxv[d]){
        maxv[d] = abs(v[e*m_dim+d]);
    }
    cout << "MAX V: "<<maxv[0]<<","<<maxv[1]<<", "<<maxv[2]<<endl;


    // (2) Update OVERALL displacements  and positions (x = x_initial + u)
    for (int i = 0; i < m_node_count * m_dim; i++) {
        u_inc[i] = dt * v[i];       // Incremental update
        //x[i] = x_initial[i] + u[i] + u_inc[i];    // Total position
        x[i] = prev_x[i] + u_inc[i];    // Total position
    }


  calcElemJAndDerivatives();
  if (!remesh_) { //Already calculated previously to account for conservation.
    CalcElemVol();  
    CalcNodalVol();
    CalcNodalMassFromVol();
  }

      
  //STRESSES CALC
  calcElemStrainRates();
  calcElemDensity();
  // if (m_dim == 3 && m_nodxelem ==4){
  //calcElemPressureANP_Nodal();
  //calcElemPressureANP();
  // }else

  if      (m_press_algorithm == 0)
    calcElemPressure();
  else if (m_press_algorithm == 1)
    calcElemPressureANP();

  calcNodalPressureFromElemental();

  
  //CalcStressStrain(dt);
  

  

  par_loop(e,m_elem_count){
    int offset = e*m_nodxelem*m_dim;  
    for (int n=0; n<m_nodxelem;n++) 
      for (int d=0;d<m_dim;d++)
        m_f_elem[offset + n*m_dim + d] = 0.0;
  }
  

    
    if (contact)
      CalcContactForces();

    bool end_it = false;
      
    solver->setZero(); //RESET K and R matrices.
    solver->beginAssembly();
    
    /////////////////////// THIS IS BEB
    //par_loop(e,m_elem_count){
    for (int e=0;e<m_elem_count;e++){
      
     // ============================================
// ELEMENT LOOP IMPLÍCITO: DEFORM + Perzyna
// ============================================

Matrix B;
B = getElemBMatrix(e); // 6 x (m_nodxelem*m_dim)
B = B * (1.0 / m_detJ[e]); // ajustar según definición B

// 1️⃣ Obtener velocidades locales del elemento
int ndof = m_nodxelem * m_dim;
Matrix vloc(ndof,1);
for (int a=0; a<m_nodxelem; a++){
    int node = getElemNode(e,a);
    for (int d=0; d<m_dim; d++){
        vloc.Set(a*m_dim + d, 0, v[node*m_dim + d]);
    }
}

// 2️⃣ Calcular tasa de deformación (Voigt)
Matrix Ddot = MatMul(B, vloc); // 6x1
double Dxx = Ddot.getVal(0,0);
double Dyy = Ddot.getVal(1,0);
double Dzz = Ddot.getVal(2,0);
double Dxy = Ddot.getVal(3,0);
double Dyz = Ddot.getVal(4,0);
double Dzx = Ddot.getVal(5,0);

// 3️⃣ Desviador y tasa equivalente
double trace = Dxx + Dyy + Dzz;
double Dev[6];
Dev[0] = Dxx - trace/3.0;
Dev[1] = Dyy - trace/3.0;
Dev[2] = Dzz - trace/3.0;
Dev[3] = Dxy;
Dev[4] = Dyz;
Dev[5] = Dzx;

double sum = Dev[0]*Dev[0] + Dev[1]*Dev[1] + Dev[2]*Dev[2]
           + 2.0*(Dev[3]*Dev[3] + Dev[4]*Dev[4] + Dev[5]*Dev[5]);
double e_dot_eq = sqrt( (2.0/3.0) * sum );
if(e_dot_eq < 1e-18) e_dot_eq = 1e-18;

// 4️⃣ Dirección del flujo plástico
Matrix n_voigt(6,1);
for(int i=0; i<6; i++){
    n_voigt.Set(i,0, Dev[i]/e_dot_eq);
}

// 5️⃣ Stress deviator y sigma_eq
Matrix s_voigt = FlatSymToVoigt(m_tau, m_dim, m_nodxelem); 
double s0=s_voigt.getVal(0,0), s1=s_voigt.getVal(1,0), s2=s_voigt.getVal(2,0);
double s3=s_voigt.getVal(3,0), s4=s_voigt.getVal(4,0), s5=s_voigt.getVal(5,0);
double sum_s = s0*s0 + s1*s1 + s2*s2 + 2.0*(s3*s3 + s4*s4 + s5*s5);
double sigma_eq = sqrt(1.5 * sum_s);

// 6️⃣ Tasa viscoplástica Perzyna
double sigma_y = mat[e]->sy0; // yield stress
double overstress = sigma_eq - sigma_y;
double dot_gamma = 0.0;
if(overstress > 0.0){
    double tau_relax = mat[e]->visc_relax_time;
    double m_perz = mat[e]->perzyna_m;
    double sigma0 = mat[e]->sy0;
    dot_gamma = (1.0 / tau_relax) * pow( overstress / sigma0, m_perz );
}

// 7️⃣ Tangente constitutiva D_gp
Matrix D_gp(6,6);
D_gp.SetZero();

// 7a) Parte desviadora estándar
for(int i=0;i<3;i++) D_gp.Set(i,i, 2.0 * dot_gamma); // normales
D_gp.Set(3,3, 2.0 * dot_gamma); // cortes
D_gp.Set(4,4, 2.0 * dot_gamma);
D_gp.Set(5,5, 2.0 * dot_gamma);

// 7b) Opcional: tangente exacta
bool useExactTangent = true; // activar/desactivar
if(useExactTangent){
    // derivada de dot_gamma respecto a sigma_eq
    double deta_de = 0.0;
    if(overstress > 0.0){
        double tau_relax = mat[e]->visc_relax_time;
        double m_perz = mat[e]->perzyna_m;
        double sigma0 = mat[e]->sy0;
        deta_de = (1.0 / tau_relax) * m_perz * pow(overstress / sigma0, m_perz - 1) / sigma0;
    }
    Matrix nT = n_voigt.getTranspose();
    Matrix outer_nn = MatMul(n_voigt, nT); // 6x6
    D_gp = D_gp + outer_nn * (2.0 * deta_de); // suma del término exacto
}

// 8️⃣ Matriz elemental K_gp
Matrix Kgp = MatMul(B.getTranspose(), MatMul(D_gp, B));
Kgp = Kgp * vol[e]; // multiplicar por volumen o peso de integración

// 9️⃣ Rhs elemental (solo -f_int)
Matrix stress_voigt = FlatSymToVoigt(m_sigma, m_dim, m_nodxelem);
Matrix fint = MatMul(B.getTranspose(), stress_voigt);
fint = fint * vol[e];

Matrix rhs = -1.0 * fint;

// 10️⃣ Ensamblaje
solver->assembleElement(e, Kgp);  // solo Kgp
solver->assembleResidual(e, rhs);         


      } // end par element loop
      
      
      //calcElemForces();
      
      if(contact)
      solver->assembleContactStiffness(1.0e8,dt);
      
       solver->finalizeAssembly();     
      
      if (contact)
        for (int n = 0; n < m_node_count*m_dim; n++)      
          solver->addToR(n,contforce[n]); //EXTERNAL FORCES



      //AFTER ASSEMBLY!
      // cout <<"K BEFORE  Dirichlet"<<endl;
      // solver->printK();
      
      //Change prescribed 
      for (int d = 0; d < m_dim; d++)
        CalcIncBCV(d);
      
      m_solver->applyDirichletBCs(); //SYMMETRY OR DISPLACEMENTS
      //cout << "Solving system"<<endl;      
      
    
