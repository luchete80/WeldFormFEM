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

using namespace std;



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

  cout << "done"<<endl;


  cout << "Imposing BCS"<<endl;
  
  
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
  
  
  double *delta_v;
  
  #ifndef BUILD_GPU
    delta_v = new double [m_dim * m_node_count];
  #else
  
  #endif
  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "------------------------MAIN LOOP---------------------------"<<endl;
  cout << "Time: "<<Time <<", End Time "<< end_t<<endl;
  while (Time < end_t) {
      
  ////// OR TIME
  if (step_count % 100 == 0){
    printf("Step %d, Time %f, End Time: %.4e, Step Time %.4e\n",step_count, Time, end_t, dt);  
    timer.click();
    //std::cout << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    std::cout << "Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    //std::cout << "CPU Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  }

  
  cout << "Storing previous values"<<endl;
  memcpy(prev_v, v, sizeof(double) * m_node_count * m_dim);
  memcpy(prev_a, a, sizeof(double) * m_node_count * m_dim);
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
  
  cout << "Imposing BCs"<<endl;
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
  cout <<"Done."<<endl;
  //cout << "----------------DISP "<<x[0]<<", "<<x[1]<<","<<x[2]<<endl;
 
  //ELEMENT PARALLEL
  
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
    cout <<"Calc derivatives and volume"<<endl;
  calcElemJAndDerivatives();
  if (!remesh_) { //Already calculated previously to account for conservation.
    CalcElemVol();  
    CalcNodalVol();
    CalcNodalMassFromVol();
  }
  cout << "Done."<<endl;
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
    cout << "Calc Strains "<<endl;
  calcElemStrainRates();
  calcElemDensity();
  // if (m_dim == 3 && m_nodxelem ==4){
  //calcElemPressureANP_Nodal();
  //calcElemPressureANP();
  // }else
    cout << "Done. Calc Pressure "<<endl;
  if      (m_press_algorithm == 0)
    calcElemPressure();
  else if (m_press_algorithm == 1)
    calcElemPressureANP();
  cout << "Done. "<<endl;
  calcNodalPressureFromElemental();

  cout << "Calc Stresses "<<endl;
  CalcStressStrain(dt);
  cout << "Done"<<endl;

  //calcArtificialViscosity(); //Added to Sigma
  
  
  // Newton-Raphson loop
  double tolerance = 1e-6;
  int max_iter = 10;
  bool converged = false;


  for (int i=0;i<m_node_count*m_dim;i++)delta_v[i]=0;
  
  cout <<"Newton Rhapson Loop"<<endl;
  ////////////////////////////////////////////////////////////
  ////////////////////////// NR LOOP /////////////////////////
  for (int iter = 0; iter < max_iter && !converged; iter++) {
  
    calcElemForces();  ///// INTERNAL FORCES
    calcElemHourglassForces();
    
    if (contact)
      CalcContactForces();


    bool end_it = false;
    
    Matrix r_global(m_nodxelem*m_dim,1);
      
    solver->beginAssembly();
    
    
    /////////////////////// THIS IS BEB
    par_loop(e,m_elem_count){
          
          // 6) Build B matrix (strain-displacement) for the element
          int tid = omp_get_thread_num();
          Matrix &B = Bmat_per_thread[tid];
          B = getElemBMatrix(e); // dimensions 6 x (m_nodxelem * m_dim)
          cout <<"Done."<<endl;
          cout << "B mat "<<endl;
          B.Print();
          cout << "m_dim "<<m_dim<<endl;
          // 7) Compute internal force: fint = V_e * B^T * σ
          cout << "Computing internal force"<<endl;
          
          Matrix stress_voigt = FlatSymToVoigt(m_sigma,m_dim,m_nodxelem);
          //CHANGE TO FORCES TO MATRIX! already calculated
          Matrix fint = MatMul(B.getTranspose(), stress_voigt); //DO NOT TRANSPOSE B DEFITELY
          /////COMPARE WITH ELEMENT FORCES
          //
          // Matrix fint;
          
          // for (int e=0;)
          // // Componentes normales
          // voigt.Set(0, 0, mat.getVal(0, 0)); // ε_xx
          // voigt.Set(1, 0, mat.getVal(1, 1)); // ε_yy
          // voigt.Set(2, 0, mat.getVal(2, 2)); // ε_zz
          // // Componentes cortantes (ingeniería: γ = 2ε)
          // voigt.Set(3, 0, mat.getVal(0, 1) + mat.getVal(1, 0)); // γ_xy = 2ε_xy
          // voigt.Set(4, 0, mat.getVal(1, 2) + mat.getVal(2, 1)); // γ_yz = 2ε_yz
          // voigt.Set(5, 0, mat.getVal(0, 2) + mat.getVal(2, 0)); // γ_xz = 2ε_xz
          // return voigt;
          
          fint = fint * vol[e];
          cout << "Calculating Kmat "<<endl;
          ///TANGENT!
          // // 8.1) Compute tangent stiffness matrix Ktan = V_e * B^T * D * B
          Matrix D(6,6);
          D =  mat[e]->getElasticMatrix();
          
          Matrix Kmat = MatMul(B.getTranspose(), MatMul(D, B));
          Kmat = Kmat * vol[e];
          // Kmat = Kmat * vol[e];
          cout << "Kmat "<<endl;
          Kmat.Print();
          
          // // 8.2) (Optional) Compute geometric stiffness Kgeo if needed


          double Ve = vol[e]; // Current volume (updated Lagrangian)


          ///////////////////////////////////////////////////
          /////////// IMPORTANT!!! --MUCH-- LESS PRODUCTS THAN: Kgeo = G^T sigma G
          // 2. Initialize Kgeo (12x12 for 4-node tetrahedron)
          //Matrix& Kgeo = *(m_Kgeo[e]);
          Matrix Kgeo(m_dim*m_node_count,m_dim*m_node_count);
          Kgeo.SetZero();

          // 3. Loop over node pairs (a, b)
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
        

          Matrix K = Kgeo + Kmat;
          
          double beta = 0.25;
          // // Add mass scaling for stability (FORGE does this)
          for (int i = 0; i < m_nodxelem*m_dim; i++) {
              //////    K[i][i] += M_diag[i] / (beta * dt * dt); // beta = 0.25 typically
              K.Set(i,i,K.getVal(i,i)+ m_mdiag[getElemNode(e,i)] / (beta * dt * dt)); // beta = 0.25 typically
          }

          Matrix R(m_dim*m_nodxelem,1);
          for (int i = 0; i < m_nodxelem * m_dim; i++) {
            //int node = getElemNode(e, i % m_nodxelem);
            R.Set(i,0,-fint.getVal(i,0)); //ADD EXTERNAL ELEMENT FORCES
          }
          ////// Residual forces (with inertial term)
          //Matrix R = f_ext - fint;
          for (int i = 0; i < m_nodxelem * m_dim; i++) {
              int node = getElemNode(e, i % m_nodxelem);
              //R[i] -= m_mdiag[node] * a[node] / (beta * dt);  // a = (v_new - v_old)/(γ*Δt)
              R.Set(i,0,R.getVal(i,0)-m_mdiag[node] * a[node] / (beta * dt));
          }
          
          solver->assembleElement(e, K);
          solver->assembleResidual(e,fint);//SHOULD BE NEGATIVE!
          
          //for (int i=0;i<m_nodxelem;i++)
            

          
          //r_global =r_global + MatMul(K,WHICH MATRIX?);

          cout << "Delta U "<<endl;
          //delta_u_e.Print();
          
      
      } // end par element loop

      m_solver->applyDirichletBCs(); //SYMMETRY OR DISPLACEMENTS
      cout << "Solving system"<<endl;      
      solver->finalizeAssembly();

      m_solver->Solve();
    
    // Update displacements and check convergence
    double max_residual = 0.0;

    for (int n = 0; n < m_node_count; n++) {
        for (int d = 0; d < m_dim; d++) {
            int idx = n * m_dim + d;
            double dv = m_solver->getU(n,d);
            delta_v[idx] += dv;
            // Track maximum residual
            max_residual = std::max(max_residual, std::abs(dv));
        }
    }
    
    // Check convergence
    if (max_residual < tolerance) {
        converged = true;
        if (step_count % 10 == 0) {
            std::cout << "NR converged in " << iter+1 << " iterations" << std::endl;
        }
    }
    
    if (iter == max_iter-1 && !converged) {
        std::cerr << "Warning: NR did not converge in " << max_iter << " iterations" << std::endl;
    }
    
  }//NR ITER 

  // After NR loop converges:
  for (int n = 0; n < m_node_count; n++) {
      for (int d = 0; d < m_dim; d++) {
          int idx = n * m_dim + d;
          // Update velocity (v_new = v_prev + Δv)
          v[idx] = prev_v[idx] +  delta_v[idx];
      }
  }

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
  /// AND AFTER REINFORCE Velocity BCs
  const double gamma = 0.5;   // Newmark parameter
  // After NR loop converges:
  for (int n = 0; n < m_node_count; n++) {
      for (int d = 0; d < m_dim; d++) {
          int idx = n * m_dim + d;

          a[idx] = (v[idx] - prev_v[idx]) / (gamma * dt) 
                 - (1.0 - gamma)/gamma * prev_a[idx];
          // Update displacement (u_new = u_prev + Δt * v_new)
          u[idx] += dt * v[idx];
          x[idx] = x[idx] + delta_v[idx]*dt;  // Ensure x is synced with u
      }
  }


  // Store values for next step
  for (int n = 0; n < m_node_count * m_dim; n++) {
      prev_v[n] = v[n];
      prev_a[n] = a[n];
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


