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
#include "tensor3.C"

using namespace std;


 tensor3 Inv(const tensor3& T) {
    tensor3 inv;
    
    // Calcular determinante
    double det = T.xx * (T.yy*T.zz - T.yz*T.zy) 
               - T.xy * (T.yx*T.zz - T.yz*T.zx) 
               + T.xz * (T.yx*T.zy - T.yy*T.zx);
    
    // Comprobar si la matriz es invertible
    if (fabs(det) < 1e-15) {
        // Matriz singular, devolver identidad como fallback
        return Identity();
    }
    
    double inv_det = 1.0 / det;
    
    // Calcular la matriz adjunta (transpuesta de la matriz de cofactores)
    inv.xx =  (T.yy*T.zz - T.zy*T.yz) * inv_det;
    inv.xy = -(T.xy*T.zz - T.zy*T.xz) * inv_det;
    inv.xz =  (T.xy*T.yz - T.yy*T.xz) * inv_det;
    
    inv.yx = -(T.yx*T.zz - T.zx*T.yz) * inv_det;
    inv.yy =  (T.xx*T.zz - T.zx*T.xz) * inv_det;
    inv.yz = -(T.xx*T.yz - T.yx*T.xz) * inv_det;
    
    inv.zx =  (T.yx*T.zy - T.zx*T.yy) * inv_det;
    inv.zy = -(T.xx*T.zy - T.zx*T.xy) * inv_det;
    inv.zz =  (T.xx*T.yy - T.yx*T.xy) * inv_det;
    
    return inv;
}
std::ostringstream m_oss;
std::string m_fname;


namespace MetFEM{

// double Domain_d::calculatePhysicsBasedTimeStep() {
    // // Base time step on physical phenomena
    // double base_dt = end_t / 100.0; // Default: 100 steps total
    
    // // Option 1: Based on maximum velocity
    // double max_vel = 0.0;
    // for(int n=0; n<m_node_count; n++) {
        // double v_mag = sqrt(v[n*m_dim+0]*v[n*m_dim+0] + 
                       // v[n*m_dim+1]*v[n*m_dim+1] + 
                       // v[n*m_dim+2]*v[n*m_dim+2]);
        // max_vel = std::max(max_vel, v_mag);
    // }
    // if(max_vel > 1e-6) {
        // double vel_based_dt = 0.1 * (domain_size / max_vel);
        // base_dt = std::min(base_dt, vel_based_dt);
    // }
    
    // // Option 2: Based on loading rate
    // if(loading_rate > 0) {
        // double loading_dt = 0.1 * (total_displacement / loading_rate);
        // base_dt = std::min(base_dt, loading_dt);
    // }
    
    // return base_dt;
// }



void host_ Domain_d::SolveImplicitGlobalMatrix(){
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


  
  
  InitValues();
  
  double3 *m_v_orig;
  if (contact){
    ////TEMP; SMOOTH LOAD(ASSUMING CONTSTANT)
    m_v_orig = new double3 [trimesh->nodecount];
    for (int n=0;n<trimesh->nodecount;n++){
      m_v_orig[n]= trimesh->node_v[n];
    }
  }
  
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
  cout << "done"<<endl;

    cout <<"VELOCITIES"<<endl;
    for (int i = 0; i < m_node_count; i++){ 
      for (int d=0;d<m_dim;d++){
      
        cout <<v[m_dim*i+d] <<", "; ///v: Current total velocity 
    }    
    cout <<endl;
    }

  
  ostringstream oss_out;
  
  
  //cout << "Calculating derivatives..."<<endl;
	#if CUDA_BUILD
	calcElemJAndDerivKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  //cout << "Calculating Volume..."<<endl;
  calcElemInitialVolKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();   
  
  calcElemDensityKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();
  
  calcElemMassMatKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();
  
  //assemblyMassMatrixKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	//cudaDeviceSynchronize();
  
  N = this->m_node_count;
	blocksPerGrid =	(N + threadsPerBlock - 1) / threadsPerBlock;
  
  CalcNodalVolKernel<<<blocksPerGrid,threadsPerBlock>>>(this);
  cudaDeviceSynchronize();
  
  CalcNodalMassFromVolKernel<<< blocksPerGrid,threadsPerBlock>>>(this);
  cudaDeviceSynchronize();
  N = this->getElemCount();
	blocksPerGrid =	(N + threadsPerBlock - 1) / threadsPerBlock;
    
  #else
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
  
  #endif
	//cout << "Done. "<<endl;

/*
  printf("INITIAL VEL\n");
  for(int e=0;e<m_elem_count;e++)
  for (int n=0;n<m_nodxelem;n++)
    printf ("elem  %d %f\n",e,getVElem(e,n,0));  
  */

  ////IMPLICIT DEFS (SAVINGS MEM)
  #ifndef BUILD_GPU
    std::vector<Matrix> Bmat_per_thread(Nproc);
    std::vector<Matrix> sig_per_thread(Nproc);
  #else
    
  #endif
  

  Time = 0.0;
  int step_count = 0;
  double tout = 0;
  
  bool remesh_ = false;
  int remesh_count = 0;
  const double RAMP_FRACTION = 1.0e-2;  // 0.1% of total time instead of 1%
  of << "t,f,area"<<endl;
  
  
  ///// SOLVER THINGS.
  Solver_Eigen *solver = new Solver_Eigen();
  m_solver = solver;
  m_solver->setDomain(this);
  m_solver->Allocate();
  
  
  double *delta_v, *x_initial;
  
  #ifndef BUILD_GPU
    delta_v   = new double [m_dim * m_node_count];
    x_initial = new double [m_dim * m_node_count];
  #else
  
  #endif
  
  double tau_old[m_elem_count*6];
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "------------------------MAIN LOOP---------------------------"<<endl;
  cout << "Time: "<<Time <<", End Time "<< end_t<<endl;
  while (Time < end_t) {

  for (int i=0;i<m_elem_count*6;i++)
    tau_old[i] = m_tau[i];
  
  ////// OR TIME
  if (step_count % 100 == 0){
    printf("Step %d, Time %f, End Time: %.4e, Step Time %.4e\n",step_count, Time, end_t, dt);  
    timer.click();
    //std::cout << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    std::cout << "Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    //std::cout << "CPU Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  }

  
  cout << "Storing previous values"<<endl;
  memcpy(prev_v,    v, sizeof(double) * m_node_count * m_dim);
  memcpy(prev_a,    a, sizeof(double) * m_node_count * m_dim);
  memcpy(x_initial, x, sizeof(double) * m_node_count * m_dim);
  cout << "Done."<<endl;
  
  cout << "Calc External Faces"<<endl;
  if (step_count % 10 == 0){
    //cout << "Calc ExtFace Areas"<<endl;
    CalcExtFaceAreas();
    //cout << "Done"<<endl;
    }
  cout << "Done"<<endl;
  
  //~ if (step_count % 50 == 0)
    //~ SearchExtNodes(); //TODO: CALCULATE ONLY AREA, NOT SEARCH AGAIN AREAS


  /////AFTER J AND DERIVATIVES
  if ( step_count % m_remesh_interval == 0 && step_count  >0 && remesh_count < m_remesh_max_count)
  //if (0) //debug
  {
    //cout << "REMAINING " <<(step_count) % m_remesh_interval<<"INTERVAL "<<m_remesh_interval<<endl;
    //cout << "step_count "<<step_count<<endl;
    double max=0.0;
    int emin;
    for (int e=0;e<m_elem_count;e++)
      if (pl_strain[e]>max){
        max = pl_strain[e];
        emin = e;
      }
      if (max>m_remesh_min_pl_strain){
  //////////////////////////// IF REMESH
      //#########################################################
      cout << "REMESHING "<< " at step "<<step_count<<endl;
      std::string ss = "in_remesh_"+std::to_string(step_count)+".vtk";
      VTKWriter writer(this, ss.c_str());
      writer.writeFile();
      
      #ifdef BUILD_REMESH
      ReMesher remesh(this);
      remesh.m_type = MMG;
      //remesh.Generate_omegah();
      remesh.Generate_mmg();
      remesh.WriteDomain(); 
      //cout << "Step "<<step_count<<endl;
      //parallel_for ()

      //TO MODIFY
      double mat_cs = sqrt(mat[0]->Elastic().BulkMod()/rho[0]);

      //cout << "Searching ext nodes "<<endl;
      SearchExtNodes(); //TODO: CALCULATE ONLY AREA, NOT SEARCH AGAIN AREAS
      //cout << "Done "<<endl;
      std::string s = "out_remesh_"+std::to_string(step_count)+".vtk";
      VTKWriter writer3(this, s.c_str());
      writer3.writeFile();
      remesh_ = true;  
      #endif
      remesh_count++;
      }
      //#########################################################
  //////////////////////////// IF REMESH

  }

  // if (!m_fixed_dt){
    // cout << "Calculate min length"<<endl;
    // double mat_cs = sqrt(mat[0]->Elastic().BulkMod()/rho[0]);
    // calcMinEdgeLength();
    // double minl = getMinLength();
    // dt = m_cfl_factor*minl/(mat_cs);
  // }
  

  //// NECESARY FOR Strain Calc
  for (int d=0;d<m_dim;d++){
    
    #ifdef CUDA_BUILD
    N = bc_count[d];
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    ImposeBCVKernel<<<blocksPerGrid,threadsPerBlock >>>(this, d);
    cudaDeviceSynchronize();
    #else
      ImposeBCV(d);
    #endif
  }

  //cout << "----------------DISP "<<x[0]<<", "<<x[1]<<","<<x[2]<<endl;
 
  //ELEMENT PARALLEL


  // Newton-Raphson loop
  double tolerance = 1e-6; //dv tol
  double ftol = 1e-6;
  int max_iter = 20;
  bool converged = false;
  double force_factor = 1.0e-3;//TO AVOID ILL CONDITIONING
  
  
  double prev_Rnorm;
  double alpha_damp= 1.0;

  ////delta_v: Pure NR correction term
  for (int i=0;i<m_node_count*m_dim;i++)delta_v[i]=0.0;
 
  
  double flat_fold[6*m_elem_count];
  
  ////////////////////////////////////////////////////////////
  ////////////////////////// NR LOOP /////////////////////////
  for (int iter = 0; iter < max_iter && !converged; iter++) {
    cout <<"ITER "<<iter<<endl;
    // cout <<"DELTA V----"<<endl;
    // for (int i = 0; i < m_node_count;i++){
      // for (int d=0;d<3;d++)
        // cout <<delta_v[m_dim*i+d]<<", ";
      // cout <<endl;
    // }
    const double gamma = 0.5;   // Newmark parameter
    // (1) Update velocities (v = prev_v + delta_v)
    for (int i = 0; i < m_node_count * m_dim; i++) {
        v[i] = prev_v[i] + delta_v[i]; ///v: Current total velocity 
    }
    
    //Rewrite BCs
    for (int d=0;d<m_dim;d++){
      
      #ifdef CUDA_BUILD
      ////REMAINS TO INIT VELOCITIES
      N = bc_count[d];
      blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
      ImposeBCVKernel<<<blocksPerGrid,threadsPerBlock >>>(this, d);
      cudaDeviceSynchronize();
      #else
      ImposeBCV(d);
      #endif
    }

    // cout <<"VELOCITIES"<<endl;
    // for (int i = 0; i < m_node_count; i++){ 
      // for (int d=0;d<m_dim;d++){
      
        // cout <<v[m_dim*i+d] <<", "; ///v: Current total velocity 
    // }    
    // cout <<endl;
    // }

    // (2) Update OVERALL displacements  and positions (x = x_initial + u)
    for (int i = 0; i < m_node_count * m_dim; i++) {
        u[i] = dt * v[i];       // Incremental update
        x[i] = x_initial[i] + u[i];    // Total position
    }

    // (3) Recompute acceleration (a) from Newmark-β
    for (int i = 0; i < m_node_count * m_dim; i++) {
        a[i] = (v[i] - prev_v[i]) / (gamma * dt) - (1.0 - gamma)/gamma * prev_a[i];
    }

    // cout <<"POSITIONS----"<<endl;
    // for (int i = 0; i < m_node_count;i++){
      // for (int d=0;d<3;d++)
        // cout <<x[m_dim*i+d]<<", ";
      // cout <<endl;
    // }

    // cout <<"DISPLACEMENTS----"<<endl;
    // for (int i = 0; i < m_node_count;i++){
      // for (int d=0;d<3;d++)
        // cout <<u[m_dim*i+d]<<", ";
      // cout <<endl;
    // }

    // cout <<"ACCELS----"<<endl;
    // for (int i = 0; i < m_node_count;i++){
      // for (int d=0;d<3;d++)
        // cout <<a[m_dim*i+d]<<", ";
      // cout <<endl;
    // }

  
  #ifdef CUDA_BUILD
  N = getElemCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;  
  
	calcElemJAndDerivKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  
	calcElemVolKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();   
  
  CalcNodalMassFromVolKernel<<< blocksPerGrid,threadsPerBlock>>>(this);
  cudaDeviceSynchronize();
  
  #else
  calcElemJAndDerivatives();
  if (!remesh_) { //Already calculated previously to account for conservation.
    CalcElemVol();  
    CalcNodalVol();
    CalcNodalMassFromVol();
  }

  #endif
   

  
  //////// END REMESH 
  ////////////////////////////////////////////
  
  #if CUDA_BUILD    
  calcElemStrainRatesKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 

  calcElemDensityKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();

  //AFTER DERIVATIVES AND RHO CALC (NEEDS JACOBIAN)
  N = getElemCount();

	// blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;  
  // calcElemMassMatKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
  // cudaDeviceSynchronize();   
  
  // //printf("CALCULATING MASS\n");
  // N = getNodeCount();
  // blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  // assemblyMassMatrixKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	// cudaDeviceSynchronize();   
 

  
    //STRESSES CALC
  N = getElemCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;  
  calcElemPressureKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
  cudaDeviceSynchronize(); 

  N = getElemCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;  
  //cout << "dt "<<dt<<endl;
  calcStressStrainKernel<<<blocksPerGrid,threadsPerBlock>>>(this, dt);
  cudaDeviceSynchronize();

  N = getNodeCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;  
  calcElemForcesKernel<<<blocksPerGrid,threadsPerBlock>>>(this);
  cudaDeviceSynchronize();

  calcElemHourglassForcesKernel<<<blocksPerGrid,threadsPerBlock>>>(this);
  cudaDeviceSynchronize();
  
  assemblyForcesKernel<<<blocksPerGrid,threadsPerBlock>>>(this);
  cudaDeviceSynchronize();
  


  #else
  //SECOND TIME
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

  
  CalcStressStrain(dt);
  



  // If not inertia terms.
  // for (int i=0;i<m_node_count;i++)
    // m_mdiag[i] = 1.0e-10;

  //calcArtificialViscosity(); //Added to Sigma
  
  
    
    calcElemForces();  ///// INTERNAL FORCES
    calcElemHourglassForces();
    
    if (contact)
      CalcContactForces();

    bool end_it = false;
    
    Matrix r_global(m_nodxelem*m_dim,1);
      
    solver->setZero(); //RESET K and R matrices.
      
    solver->beginAssembly();
    
    /////////////////////// THIS IS BEB
    par_loop(e,m_elem_count){
          cout << "Element "<<e<<endl;
          
          // 6) Build B matrix (strain-displacement) for the element
          int tid = omp_get_thread_num();

          Matrix &B = Bmat_per_thread[tid];
          //// HERE B is in fact BxdetJ
          B = getElemBMatrix(e); // dimensions 6 x (m_nodxelem * m_dim)
          B = B *(1.0/m_detJ[e]);

          //cout << "B mat "<<endl;
          //B.Print();

          // 7) Compute internal force: fint = V_e * B^T * σ

          Matrix stress_voigt = FlatSymToVoigt(m_sigma,m_dim,m_nodxelem);
          //CHANGE TO FORCES TO MATRIX! already calculated
          Matrix fint = MatMul(B.getTranspose(), stress_voigt); //DO NOT TRANSPOSE B DEFITELY
 
          fint = fint * vol[e];

          ///TANGENT!
          // // 8.1) Compute tangent stiffness matrix Ktan = V_e * B^T * D * B
          Matrix D(6,6);
          D =  mat[e]->getElasticMatrix();
          
          Matrix Kmat = MatMul(B.getTranspose(), MatMul(D, B));
          Kmat = Kmat * (1.0/6.0*m_detJ[e]); // B is B x detJ

          double Ve = vol[e]; // Current volume (updated Lagrangian)

          ///////////////////////////////////////////////////
          /////////// IMPORTANT!!! --A LOT-- FASTER (LESS PRODUCTS) THAN: Kgeo = G^T sigma G
          // 2. Initialize Kgeo (12x12 for 4-node tetrahedron)
          //Matrix& Kgeo = *(m_Kgeo[e]);
          Matrix Kgeo(m_dim*m_nodxelem,m_dim*m_nodxelem);
          Kgeo.SetZero();
          
          // // 3. Loop over node pairs (a, b)
          // // REMEMBER DERIVATIVES ARE AFFECTED BY DETJ
          for (int a = 0; a < 4; ++a) {
            // ∇Nᵃ in current config (∂Nᵃ/∂x, ∂Nᵃ/∂y, ∂Nᵃ/∂z)
            Matrix grad_a(3, 1);
            grad_a.Set(0, 0, getDerivative(e, 0, 0, a)); // ∂N/∂x
            grad_a.Set(1, 0, getDerivative(e, 0, 1, a)); // ∂N/∂y
            grad_a.Set(2, 0, getDerivative(e, 0, 2, a)); // ∂N/∂z

            for (int b = 0; b < 4; ++b) {
              // ∇Nᵇ in current config
              Matrix grad_b(3, 1);
              grad_b.Set(0, 0, getDerivative(e, 0, 0, b));
              grad_b.Set(1, 0, getDerivative(e, 0, 1, b));
              grad_b.Set(2, 0, getDerivative(e, 0, 2, b));

              // Compute K_geo(a,b) = (∇Nᵃ)ᵀ · σ · ∇Nᵇ * Ve
              Matrix sigma_grad_b = MatMul(FlatSymToMatrix(m_sigma), grad_b); // σ · ∇Nᵇ (3x1)
              Matrix kab = MatMul(grad_a.getTranspose(), sigma_grad_b); // 1x1 scalar
              double k_ab = kab.getVal(0, 0) * Ve;

              // Fill 3x3 block (assumes 3 DOF per node)
              for (int i = 0; i < 3; ++i) {
                Kgeo.Set(3*a + i, 3*b + i, Kgeo.getVal(3*a + i, 3*b + i) + k_ab);
              }
            }
          }

          cout << "Kgeo row col "<<Kgeo.m_row<<","<<Kgeo.m_col<<endl;
          Kgeo = Kgeo * (1.0/(6.0*m_detJ[e]));
          cout <<"sum "<<endl;
          cout << "Kmat row col "<<Kmat.m_row<<","<<Kmat.m_col<<endl;
          Matrix K = Kgeo + Kmat;
          cout << "done"<<endl;
          K = K*dt;
          
          double beta = 0.25;
          // // Add mass scaling for stability (FORGE does this)
          for (int i = 0; i < m_nodxelem; i++) {  // Loop over element nodes
              int node = getElemNode(e, i);        // Get global node number
              for (int d = 0; d < m_dim; d++) {   // Loop over dimensions (x,y,z)
                  int idx = i*m_dim + d;           // Local DOF index
                  double mass_term = m_mdiag[node] / (beta * dt);  //kg/s = (kgxm/s2) x s/m = N/m x s
                  K.Set(idx, idx, (K.getVal(idx, idx) + mass_term) *(1.0 + 1.0e-8) ); //ALSO ADDED DIAG REGULARIZATION
              }
          }
          //cout <<"CHECKING INTERNAL FORCES"<<endl;

          Matrix R(m_dim*m_nodxelem,1);
          for (int i = 0; i < m_nodxelem; i++) {
            //int node = getElemNode(e, i % m_nodxelem);
            for (int d=0;d<m_dim;d++){
            //cout << "NODE, DIM "<<i<<","<<d<<", fint mat"<<fint.getVal(m_dim*i+d,0)<<", fel "<<m_f_elem[i*m_dim+d]<<endl;
            //R.Set(i,0,-fint.getVal(m_dim*i+d,0)); //ADD EXTERNAL ELEMENT FORCES
            R.Set(m_dim*i+d,0,-m_f_elem[i*m_dim+d]/*+m_f_elem_hg [offset + i*m_dim + d]*/); //ADD EXTERNAL ELEMENT FORCES
            }
          }

          // cout <<"INTERTIA TERMS OF RESIDUAL"<<endl;
          // ////// Residual forces (with inertial term)
          // //Matrix R = f_ext - fint;
          // for (int i = 0; i < m_nodxelem; i++) {
              // int node = getElemNode(e, i);
              // for (int d=0;d<m_dim;d++){
              // //R[i] -= m_mdiag[node] * a[node] / (beta * dt);  // a = (v_new - v_old)/(γ*Δt)
                // cout << "Node DIM "<<node<<","<<d<<", "<<"R Orig"<<R.getVal(gdof,0)<<"Inertia"<<-m_mdiag[node] * a[gdof]<<endl;
                // R.Set(i,0,R.getVal(gdof,0)-m_mdiag[node] * a[gdof]);

              // }
          // }
          
          solver->assembleElement(e, K);
          solver->assembleResidual(e,R);//SHOULD BE NEGATIVE!  
      
        //cout << "Element R "<<endl;
        //R.Print();
      } // end par element loop
      
      

      for (int n = 0; n < m_node_count*m_dim; n++)      
        solver->addToR(n,contforce[n]); //EXTERNAL FORCES


      for (int n = 0; n < m_node_count; n++){   
        for (int d=0;d<m_dim;d++){        
          int gdof = m_dim*n+d;
          solver->addToR(gdof,-m_mdiag[n] * a[gdof] ); //EXTERNAL FORCES
        }
      }    
      // cout <<"R AFTER INERTIA AND CONTACT"<<endl;
      // solver->printR();
      
      solver->finalizeAssembly();
      //AFTER ASSEMBLY!
      // cout <<"K BEFORE  Dirichlet"<<endl;
      // solver->printK();

      
      m_solver->applyDirichletBCs(); //SYMMETRY OR DISPLACEMENTS
      //cout << "Solving system"<<endl;      
      
      m_solver->Solve();
    
    // Update displacements and check convergence
    double max_residual = 0.0;
    double max_f_residual = 0.0;


    // if (iter >0){
      // double residual_ratio = m_solver->getRNorm() / prev_Rnorm; //Ratio is larger than 1 if is performing bad
      // alpha_damp = std::min(1.0, 1.0 / (1.0 + 1000.0 * residual_ratio));    
      // cout << "Using alpha damp: "<<alpha_damp<<endl;
    // }
    prev_Rnorm = m_solver->getRNorm();    
    
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
    
      
    cout << "MAX Residuals, DV: "<< max_residual<<
    ", DF ABS: "<<max_f_residual<<endl;
    
    cout << "Urel "<<max_residual/m_solver->getUNorm()<<"U Norm "<<m_solver->getUNorm()<<endl;
    cout << "frel "<<max_f_residual/m_solver->getRNorm()<<"R Norm "<<m_solver->getRNorm()<<endl;
    
    // Check convergence
    if (max_residual < tolerance && max_f_residual < ftol) {
        converged = true;
        if (step_count % 10 == 0) {
            std::cout << "NR converged in " << iter+1 << " iterations" << std::endl;
        }
    }
    
    if (iter == max_iter-1 && !converged) {
        std::cerr << "Warning: NR did not converge in " << max_iter << " iterations" << std::endl;
    }
    
    //SWITCH BACK TO PREVIOUS STRESS STATE, WHICH IS TAKEN BY calcStressStrain
    if (!converged){
      for (int i=0;i<m_elem_count*6;i++)
        m_tau[i] = tau_old[i];
    }

    
  }//NR ITER 
  /////////////////////////////////////////////////////////////////////////////////////////////////

  // ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
  // //ToFlatSymPtr(Strain, 	strain,6*i);		
  // ToFlatSymPtr(ShearStress, m_tau, offset_t);
      
  // After NR converges:
  for (int i = 0; i < m_node_count * m_dim; i++) {
      prev_v[i] = v[i];  // Save converged velocity
      prev_a[i] = a[i];  // Save acceleration
  }
  
  ///assemblyForces(); 
  //ApplyGlobalSprings();

  
  #endif
  
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
    
  
  if (contact){
    //MeshUpdate(this->trimesh,dt);
  #ifdef CUDA_BUILD  
  #else
    trimesh->Move( dt);
    trimesh->CalcCentroids();
    trimesh->CalcNormals();

    trimesh->UpdatePlaneCoeff();
  #endif
  }
  if (remesh_){
   //printf("DISPLACEMENTS\n");
   //printVec(this->u);       
       std::string s = "out_remesh_after1_"+std::to_string(step_count)+".vtk";
      VTKWriter writer3(this, s.c_str());
      writer3.writeFile();   
    }

  ///// AFTER CONTACT (FOR THE cont_cond)
  if(m_thermal){
    //calcInelasticHeatFraction(); //Before thermal calc
    ThermalCalcs(); //m_dTedt[e1n1 e1n2 e1n3 e1n4 _ e2n1 ..]

}

 
  if (Time>=tout){
    string outfname = "out_" + std::to_string(Time) + ".vtk";
    timer.click();

    ostringstream oss_out;
    oss_out << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    oss_out << "CPU Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    oss_out << "Plastic Strain energy "<<m_pl_energy<<endl;

    of <<std::scientific<<std::setprecision(6)<< Time ;

    if (contact){
    calcContactForceFromPressure();
    
    of <<std::scientific<<std::setprecision(6)<<", "<<
                                                    //trimesh->react_p_force[m]<<
                                                    trimesh->react_p_force[0]<<", "<<
                                                    //trimesh->react_force[0].z<<","<<
                                                    trimesh->cont_area;
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
          
          

    
  
  Time += dt;
  step_count++;
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
    /*
  calcElemStrainRates();
  
   printf("DISPLACEMENTS\n");
   printVec(this->u);   

  printf("VELOCITIES\n");
  printVec( this->v);

  printf("ACCEL\n");
	printVec(this->a); 
  
  printf("FORCES\n");
  printVec(this->m_fi);
*/
  // printf("STRESSES\n");
  // printSymmTens(this->m_sigma);

  // printf("SHEAR STRESS\n");
  // printSymmTens(this->m_tau);

  // printf("STRAIN RATES\n");
  // printSymmTens(this->m_str_rate);
  
  // printf("ROT RATES\n");
  // printSymmTens(this->m_rot_rate);
  
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


