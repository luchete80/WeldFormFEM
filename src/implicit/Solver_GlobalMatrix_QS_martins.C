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
      
    // // // // Dentro de loop elementos
    // // // Matrix B = getElemBMatrix(e); B = B*(1.0/m_detJ[e]);
    // // // Matrix vloc = extract_vloc(e);
    // // // Matrix Ddot = MatMul(B, vloc); // 6x1
    // // // compute Dev, e_dot_eq, n_voigt from Ddot
    // // // compute sigma_y = ... // yield
    // // // compute sigma_eq_trial from s_trial if using trial, else use last sigma

    // // // compute dot_gamma (Perzyna) OR check f_trial
    // // // if elastic:
       // // // D_gp = mat[e]->getElasticMatrix();
    // // // else if radial-return:
       // // // D_gp = getConsistentPlasticTangentMatrix(...)
    // // // else if perzyna-visco:
       // // // build C_visc as shown (Pdev + outer_nn term)

    // // // Matrix Kgp = MatMul(B.getTranspose(), MatMul(D_gp, B)) * vol[e];
    // // // Matrix b_gp = MatMul(Kgp, vloc);
    // // // Pstar_elem += b_gp * (1.0/e_dot_eq) * vol[e];

    // // // Matrix bbT = MatMul(b_gp, b_gp.getTranspose());
    // // // double dsig_dedot = ... // model-specific
    // // // double coef = ( dsig_dedot*(1.0/e_dot_eq) - sigma_bar/(e_dot_eq*e_dot_eq) )*(1.0/e_dot_eq);
    // // // H_elem += bbT * (coef * vol[e]);

    // // // Matrix Bvol = MatMul(Cmat, B); // 1xndof
    // // // Q_elem += MatMul( Bvol.getTranspose(), Bvol ) * vol[e];

    // // // // Form Aelem and rhs, assemble



          //Matrix &B = Bmat_per_thread[tid];
          Matrix B;
          int offset_t = e*6;
          
          //// HERE B is in fact BxdetJ
          B = getElemBMatrix(e); // dimensions 6 x (m_nodxelem * m_dim)
          B = B *(1.0/m_detJ[e]);

          Matrix Ddot = FlatSymToVoigt(m_str_rate, e*6, m_dim, m_nodxelem);

          int ndof = m_nodxelem * m_dim;
          Matrix vloc(ndof, 1);
          Matrix delta_vloc(ndof, 1);
          for (int a=0;a<m_nodxelem;a++){
            int node = getElemNode(e,a);
            for (int d=0; d<m_dim; d++){
              vloc.Set(a*m_dim + d, 0, v[node*m_dim + d]);
              delta_vloc.Set(a*m_dim + d, 0, delta_v[node*m_dim + d]);
            }
          }          
                    //Matrix Ddot = MatMul(B, vloc); // 6x1
          //Strain rate
          /////D - Ddot_dev y e_dot_eq (tasa equivalente)
          //Matrix Ddot = MatMul(B, vloc); // 6x1
          double Dxx = Ddot.getVal(0,0);
          double Dyy = Ddot.getVal(1,0);
          double Dzz = Ddot.getVal(2,0);
          double Dxy = Ddot.getVal(3,0); // engineering gamma_xy
          double Dyz = Ddot.getVal(4,0);
          double Dzx = Ddot.getVal(5,0);

          // volumetric rate
          double trace = Dxx + Dyy + Dzz;
          double Dvol = trace; // since B gives engineering strain rates (norms consistent)

          // desviador (Voigt ingenieril): normales - 1/3 trace, cortes dejar igual
          double Dev[6];
          Dev[0] = Dxx - trace/3.0;
          Dev[1] = Dyy - trace/3.0;
          Dev[2] = Dzz - trace/3.0;
          Dev[3] = Dxy;
          Dev[4] = Dyz;
          Dev[5] = Dzx;

          // norma equivalente: e_dot_eq = sqrt(2/3 * (Dev : Dev_tensor))
          double sum = Dev[0]*Dev[0] + Dev[1]*Dev[1] + Dev[2]*Dev[2]
                     + 2.0*(Dev[3]*Dev[3] + Dev[4]*Dev[4] + Dev[5]*Dev[5]); // engineering
          double e_dot_eq = sqrt( (2.0/3.0) * sum );
          if (e_dot_eq < 1e-18) e_dot_eq = 1e-18;
          
          ///E) n_COIGT FLUX DIRECTION
          Matrix n_voigt(6,1);
          n_voigt.Set(0,0, Dev[0]/e_dot_eq);
          n_voigt.Set(1,0, Dev[1]/e_dot_eq);
          n_voigt.Set(2,0, Dev[2]/e_dot_eq);
          // para ingenieril, las componentes de corte deben ser multiplicadas por 1 si Dev ya contiene gamma; 
          // si usas convención de n_t = (3/2) s / sigma_eq en retorno radial, ajusta como hiciste
          n_voigt.Set(3,0, Dev[3]/e_dot_eq);
          n_voigt.Set(4,0, Dev[4]/e_dot_eq);
          n_voigt.Set(5,0, Dev[5]/e_dot_eq);

          
          //cout << "B mat "<<endl;
          //B.Print();


          /////F) Sigma Eq
          Matrix s_voigt = FlatSymToVoigt(m_tau, e*6,m_dim, m_nodxelem); // o usa offset
          // implementar norma: sigma_eq = sqrt(3/2 * s_dev : s_dev)
          double s0 = s_voigt.getVal(0,0);
          double s1 = s_voigt.getVal(1,0);
          double s2 = s_voigt.getVal(2,0);
          double s3 = s_voigt.getVal(3,0); // tau_xy (engineering)
          double s4 = s_voigt.getVal(4,0);
          double s5 = s_voigt.getVal(5,0);
          double sum_s = s0*s0 + s1*s1 + s2*s2 + 2.0*(s3*s3 + s4*s4 + s5*s5);
          double sigma_eq = sqrt(1.5 * sum_s);
          
          tensor3 s_test = FromFlatSym(m_tau, offset_t );
          
          //sigma_eq = K * pow(e_dot_eq, m);
          
          ///////G) Perzyna: eta_eff o dot_gamma (tasa viscoplástica) y d(eta)/de si hace falta para tangente
          
          double sigma_y = CalcHollomonYieldStress(pl_strain[e], mat[e]);/* tu CalcHollomonYieldStress(pl_strain[e], mat[e]) */;

          double overstress = sigma_eq - sigma_y;
          double dot_gamma = 0.0;
          double tau_relax = mat[e]->visc_relax_time; // define en material (ej.)
          double m_perz = mat[e]->perzyna_m;
          
          if (overstress > 0.0) {
            double sigma0 = mat[e]->sy0; // escala
            dot_gamma = (1.0 / tau_relax) * pow( overstress / sigma0, m_perz );
            
            double dep = dt * dot_gamma;
            pl_strain[e] += dep; // dep calculado según Perzyna
            
            
          } else {
            dot_gamma = 0.0;
          }

          double eta = 0.5 * sigma_eq / e_dot_eq; // careful, only meaningful if sigma_eq defined from current state
          
          //tensor3 Sigma = shear_stress;
          
          ////// RIGID PLASTIC
          ////s=2ηε˙dev​,η=2ε˙eq​σeq​​
          s_voigt.Set(0,0, 2.0*eta*Dev[0]);
          s_voigt.Set(1,0, 2.0*eta*Dev[1]);
          s_voigt.Set(2,0, 2.0*eta*Dev[2]);
          s_voigt.Set(3,0, 2.0*eta*Dev[3]);
          s_voigt.Set(4,0, 2.0*eta*Dev[4]);
          s_voigt.Set(5,0, 2.0*eta*Dev[5]);
          ToFlatMat(s_voigt,m_tau,e*6);

          //UPDATE
          sigma_eq = sigma_y * pow(tau_relax * e_dot_eq, 1.0/m_perz);
          
          double K  = mat[e]->Elastic().BulkMod();
          p[e] = K * Dvol;
          
          //ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM	
          //ToFlatSymPtr(shear_stress, m_tau, offset_t); 
          //ToFlatSymPtr(Strain_pl_incr, m_strain_pl_incr, offset_t);
   
          /////////////////////////////////////////////////////////////
          // H) D_gp — Tangente viscoplástica DIAGONAL (Martins)
          /////////////////////////////////////////////////////////////

          Matrix D_gp(6,6);
          D_gp.SetZero();

          // normales
          D_gp.Set(0,0, 2.0/3.0);
          D_gp.Set(1,1, 2.0/3.0);
          D_gp.Set(2,2, 2.0/3.0);

          // cortes (ingenieril)
          D_gp.Set(3,3, 1.0/3.0);
          D_gp.Set(4,4, 1.0/3.0);
          D_gp.Set(5,5, 1.0/3.0);

          // factor común
          D_gp = D_gp * (2.0 * eta);

          /////////////////////////////////////////////////////////////
          // I) Kgp = Bᵀ · D_gp · B · vol
          /////////////////////////////////////////////////////////////

          Matrix Kgp = MatMul(B.getTranspose(), MatMul(D_gp, B));
          Kgp = Kgp * vol[e];

          /////////////////////////////////////////////////////////
          /// J) b_gp = Kgp * vloc (ndof×1)

          Matrix b_gp = MatMul(Kgp, vloc);
          
          /////// K) Pstar (vector ndof×1)
          
          Matrix Pstar(ndof,1); Pstar.SetZero();
          Pstar = Pstar + b_gp * ( (1.0 / e_dot_eq) * vol[e] );

          //////////////// L) H_elem (ndof×ndof) — integrar coef * (b b^T)
          //////////////// Coef según Martins / tu imagen (simplificada operativa):

          // compute sigma_bar (||s||) from s_dev (use shear stress trial or current)

          // derivative dsigma/de (depends on model)  — for Norton/Kelvin you can get analytic
          //HOLLOMON CASE 
          double dsig_dedot = /* analytic or approximate */ 0.0;

          // coef
          double coef = ( dsig_dedot * (1.0 / e_dot_eq) - sigma_eq / (e_dot_eq*e_dot_eq) ) * (1.0 / e_dot_eq);

          // H contribution:
          Matrix H_elem = coef * MatMul( b_gp, b_gp.getTranspose() ); // ndof x ndof
        
          
          ///////////////M) Q_elem (ndof×ndof) — volumetric coupling
          Matrix Cmat(1,6); Cmat.SetZero();
          Cmat.Set(0,0,1.0); Cmat.Set(0,1,1.0); Cmat.Set(0,2,1.0);

          Matrix Bvol = MatMul(Cmat, B); // 1 x ndof
          Matrix Qgp = MatMul( Bvol.getTranspose(), Bvol ); // ndof x ndof
          Qgp = Qgp * vol[e];

          //Q_elem = Q_elem + Qgp;
          Matrix Q_elem =  Qgp;
          
          //// N) Assemble elemental matrix A_elem = H_elem + Kgp + Kgp * Q_elem ???
          Matrix Aelem = Kgp + H_elem + Q_elem;    

          // Rhs elemental: R = f_ext_elem - f_int_elem - sigma_bar*Pstar - Kgp * (Q_elem * vloc)
          Matrix stress_voigt = FlatSymToVoigt(m_sigma,e*6,m_dim,m_nodxelem);
          //CHANGE TO FORCES TO MATRIX! already calculated
          Matrix fint = MatMul(B.getTranspose(), stress_voigt); //DO NOT TRANSPOSE B DEFITELY
          fint = fint * vol[e];
          Matrix rhs = /*f_ext_elem */-1.0* fint; // f_ext_elem may be zeros except contact
          rhs = rhs - Pstar * sigma_eq; // Pstar scalar*vector
          Matrix temp = MatMul( Q_elem, vloc );
          rhs = rhs - MatMul(Kgp, temp);
          
          /// ASSEMBLY 
          solver->assembleElement(e, Aelem);
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
      
      m_solver->Solve();
    
    // Update displacements and check convergence
    double max_residual = 0.0;
    double max_f_residual = 0.0;
    
    double max = 0.0;
    for (int e=0;e<m_elem_count;e++)
      if (pl_strain[e]>max){
        max = pl_strain[e];
      }
     
    cout << "Max plastic strain: "<<max<<endl;

    // if (iter >0){
      // double residual_ratio = m_solver->getRNorm() / prev_Rnorm; //Ratio is larger than 1 if is performing bad
      // alpha_damp = std::min(1.0, 1.0 / (1.0 + 1000.0 * residual_ratio));    
      // cout << "Using alpha damp: "<<alpha_damp<<endl;
    // }
    prev_Rnorm = m_solver->getRNorm();    

    alpha_damp = 1.0;
    if (iter > 0) {
        double residual_ratio = m_solver->getRNorm() / prev_Rnorm;
        
        if (residual_ratio > 2.0) {
            // Divergencia - reducir paso
            alpha_damp = 0.5;
        } else if (residual_ratio > 1.2) {
            // Convergencia lenta
            alpha_damp = 0.8;
        } else if (residual_ratio < 0.8) {
            // Buena convergencia - aumentar paso
            alpha_damp = std::min(1.5, alpha_damp * 1.2);
        } else {
            // Convergencia estable
            alpha_damp = 1.0;
        }
        
        // Límites más amplios
        alpha_damp = std::max(0.3, std::min(2.0, alpha_damp));
    }
    
    for (int n = 0; n < m_node_count; n++) {
        for (int d = 0; d < m_dim; d++) {
            int idx = n * m_dim + d;
            double dv = m_solver->getU(n,d);
            double df = m_solver->getR(n,d);
            delta_v[idx] += alpha_damp*dv;
            // Track maximum residual
            max_residual = std::max(max_residual, std::abs(dv));
            max_f_residual = std::max(max_residual, std::abs(df));
        }
    }
      
      double F_ref = 0.0;
      // 1. Fuerzas externas aplicadas (contacto, cargas, etc.)
      for (int n = 0; n < m_node_count * m_dim; n++) {
          F_ref = std::max(F_ref, std::abs(contforce[n])); // Fuerzas de contacto
          // Agregar otras fuerzas externas si las tienes
      }
      
      // 2. Fuerzas internas representativas (opcional pero recomendado)
      for (int e = 0; e < m_elem_count; e++) {
          for (int i = 0; i < m_nodxelem * m_dim; i++) {
              F_ref = std::max(F_ref, std::abs(m_f_elem[e * m_nodxelem * m_dim + i]));
          }
      }
      
    
      
    cout << "MAX Residuals, DV: "<< max_residual<<
    ", DF ABS: "<<max_f_residual<<", DF rel "<< max_f_residual/F_ref<<endl;
    
    // cout << "Urel "<<max_residual/m_solver->getUNorm()<<"U Norm "<<m_solver->getUNorm()<<endl;
    // cout << "frel "<<max_f_residual/m_solver->getRNorm()<<"R Norm "<<m_solver->getRNorm()<<endl;
    
    // Check convergence
    if (max_residual < tolerance && max_f_residual/F_ref < ftol) {
        if (iter>0)
          converged = true;
        end = true;
        //if (step_count % 10 == 0) {
            std::cout << "NR converged in " << iter+1 << " iterations" << std::endl;
        //}
    }
    
    if(max_residual>1.0e3 || iter >max_iter){
      
      end= true;
      dt /=2.0;
      converged = false;
      conv_iter = 0;
    }
    
    //if (iter >0) converged = true;
    
    if (iter == max_iter-1 && !converged) {
        std::cerr << "Warning: NR did not converge in " << max_iter << " iterations" << std::endl;


    }
    
    //~ //SWITCH BACK TO PREVIOUS STRESS STATE, WHICH IS TAKEN BY calcStressStrain
    if (!converged){
      //~ for (int i=0;i<m_elem_count*6;i++)
        //~ m_tau[i] = tau_old[i];
      memcpy(pl_strain, pls_old, sizeof(double) * m_elem_count);
      memcpy(m_sigma, sig_old, sizeof(double) * m_elem_count * m_gp_count * 6);
      memcpy(m_tau, tau_old, sizeof(double) * m_elem_count * m_gp_count * 6);
    
    }
    
    else {
      
          conv_iter++;
      }
  
    iter++;

    
    if (converged && conv_iter > 2){
        dt *=2.0;
    }
    

  }//NR ITER //////////////////////////////// INNER LOOP
  

  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////

  // ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
  // //ToFlatSymPtr(Strain, 	strain,6*i);		
  // ToFlatSymPtr(ShearStress, m_tau, offset_t);
    
    
    
  if (converged){
    // After NR converges:
    for (int i = 0; i < m_node_count * m_dim; i++) {
        prev_x[i] = x[i];  // Save converged velocity
        prev_v[i] = v[i];  // Save converged velocity
        
        u[i] += u_inc[i];
    }
    
      calcElemJAndDerivatives();
  if (!remesh_) { //Already calculated previously to account for conservation.
    CalcElemVol();  
    CalcNodalVol();
    CalcNodalMassFromVol();
  }
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

    
    CalcStressStrain(dt);

    N = getNodeCount();
    //printf("Correction\n");	
    #ifdef CUDA_BUILD
    
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    
    #else
    
    if (contact){
      //if (Time > RAMP_FRACTION*end_t)
      //ApplyGlobalDamping(0.02);
    }
    
    if (remesh_){
      //if (Time > RAMP_FRACTION*end_t)
      ApplyGlobalDamping(m_remesh_damp_vel);
    }

    //ApplyGlobalDamping(0.1);
    #endif



    // if (contact){
      // double f =1.0;

      // if(Time < RAMP_FRACTION*end_t) {
          // f = pow(Time/(RAMP_FRACTION*end_t), 0.5);  // Square root for smoother start
      // } else {
          // f = 1.0;
      // }
      // for (int n=0;n<trimesh->nodecount;n++){
        // trimesh->node_v[n] = f*m_v_orig[n];
      // }
        // //cout << "Node 0 v"<<(trimesh->node_v[0]).z<<endl;
    // }

    
  
    string outfname = "out_" + std::to_string(iter) + ".vtk";
    timer.click();

    ostringstream oss_out;
    oss_out << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    oss_out << "CPU Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    oss_out << "Plastic Strain energy "<<m_pl_energy<<endl;

    of <<std::scientific<<std::setprecision(6)<< iter ;

    if (contact){
    calcContactForceFromPressure();
    
    of <<std::scientific<<std::setprecision(6)<<", "<<
                                                    //trimesh->react_p_force[m]<<
                                                    trimesh->react_p_force[0]<<", "<<
                                                    //trimesh->react_force[0].z<<","<<
                                                    trimesh->cont_area;
                                                    


    
  } else{ ///// NOT CONVERGED
    for (int i = 0; i < m_node_count * m_dim; i++) {
        x[i] = prev_x[i];  // Save converged velocity
        v[i] = prev_v[i];  // Save converged velocity

    }
       
    
    
    }//Not converged
  ///assemblyForces(); 
  //ApplyGlobalSprings();

    
 
  if (Time>=tout){
    // string outfname = "out_" + std::to_string(Time) + ".vtk";
    // timer.click();

    // ostringstream oss_out;
    // oss_out << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    // oss_out << "CPU Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    // oss_out << "Plastic Strain energy "<<m_pl_energy<<endl;

    // of <<std::scientific<<std::setprecision(6)<< Time ;

    // if (contact){
    // calcContactForceFromPressure();
    
    // of <<std::scientific<<std::setprecision(6)<<", "<<
                                                    // //trimesh->react_p_force[m]<<
                                                    // trimesh->react_p_force[0]<<", "<<
                                                    // //trimesh->react_force[0].z<<","<<
                                                    // trimesh->cont_area;
    } else{
   bool is_elem_sum[m_elem_count];

   //double pxa_el[m_elem_count];
  double zmax = 0.0;
      for (int i=0;i<m_node_count;i++)
        if (getNodePos3(i).z>zmax)
          zmax = getNodePos3(i).z;
        
    int ecount = 0;
   double area = 0.0;
   
         for (int e=0;e<m_elem_count;e++)is_elem_sum[e]=false;
      double cfsum = 0.0;
      for (int i=0;i<m_node_count;i++){
        if ((getNodePos3(i).z-zmax)*(getNodePos3(i).z-zmax)<1.0e-5){
          for (int ne=0; ne<m_nodel_count[i];ne++) {
            int e   = m_nodel     [m_nodel_offset[i]+ne]; //Element
            if (!is_elem_sum[e]){
              //pxa_el[e]+=p[e]*m_elem_area[e];
              is_elem_sum[e]=true;
              cfsum += p[e]*m_elem_area[e];
              area+=m_elem_area[e];
              ecount++;
            }
              //~ if (!is_node_sum[i]){
                //~ area+=node_area[i];
                //~ }
          }//nodel
        }//mesh in contact
        
      }//Node count    
      cout << "Cont Elements"<<ecount<<endl;
      of <<std::scientific<<std::setprecision(6)<<", "<<cfsum;
    } //NOT CONTACT, TO DELETE
    
    
  double max[]={0.0,0.0,0.0};

     for (int e=0;e<m_node_count;e++)
        for (int d=0;d<3;d++)
              if (this->u[e*m_dim+d]*this->u[e*m_dim+d]>max[d]){
          max[d] = this->u[e*m_dim+d];

      } 
  
  cout << "MAX DISP "<<sqrt(max[0])<< " "<<sqrt(max[1])<< " "<<sqrt(max[2])<<endl;                                                     

    
    cout << oss_out.str();
    if (!out_file.is_open()) {
        std::cerr << " out_file is not open! Cannot write." << std::endl;
    } else {
        out_file << oss_out.str();
        out_file.flush();
    }
    
    
    //}
    of <<endl;
    #ifndef CUDA_BUILD
    VTKWriter writer2(this, outfname.c_str());
    writer2.writeFile();
    #endif
    tout +=m_dtout;
  }
      ///////DEBGUG
      //std::string s = "out_step_"+std::to_string(step_count)+".vtk";
      //VTKWriter writer3(this, s.c_str());
      //writer3.writeFile();
          
          

    
  if (converged){
    Time += dt;
    step_count++;
  
  }
  remesh_ = false;    
  }// WHILE LOOP


  //////////////////////////// IF REMESH
  #ifdef BUILD_REMESH
  ReMesher remesh(this);
  
  remesh.Generate_mmg();
  remesh.m_type = MMG;
  #endif
  //////////////////////////////////////
  

  #ifdef CUDA_BUILD
  cudaMemcpy(x_h, x, 3*sizeof(double) * m_node_count, cudaMemcpyDeviceToHost);		
  


  #else
  double max[]={0.0,0.0,0.0};
  
     for (int e=0;e<m_node_count;e++)
        for (int d=0;d<3;d++)
              if (this->u[e*m_dim+d]>max[d]){
          max[d] = this->u[e*m_dim+d];

      } 
  
  cout << "MAX DISP "<<max[0]<< " "<<max[1]<< " "<<max[2]<<endl;

  #endif
  cout << "Writing output "<<endl;
  //VTUWriter writer(this, "out.vtu");
  //writer.writeFile();

  
  #ifndef CUDA_BUILD
  //cout << "Writing output"<<endl;
  VTKWriter writer2(this, "out.vtk");
  writer2.writeFile();
  #endif
  cout << "Done."<<endl;

  ///// IF REMESH
  /////#####################
  //remesh.WriteDomain();
  //calcElemJAndDerivatives();    
  
  //////////////////////////////////////
  
  VTKWriter writer3(this, "out_remesh.vtk");
  writer3.writeFile();
  
  //AFTER WRITE

  timer.stop();
  std::cout << "Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  of.close();
  }//SOLVE
    
};

