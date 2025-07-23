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


// Matrix ComputeGlobalResidual() {
    // Matrix r_global(m_node_count * m_dim, 1);
    // r_global.SetZero();

    // par_loop(e, m_elem_count) {
        // // 1. Compute element contribution
        // Matrix B = getElemBMatrix(e);
        // Matrix f_int_e = MatMul(B.Transpose(), stress_voigt) * vol[e];
        
        // // 2. Compute K_e * u_e (stiffness action)
        // Matrix u_e = gatherElementDisplacements(e);
        // Matrix K_e = MatMul(B.Transpose(), MatMul(D, B)) * vol[e];
        // Matrix Ku_e = MatMul(K_e, u_e);
        
        // // 3. Scatter to global residual
        // scatterAdd(r_global, Ku_e - f_int_e, e);
    // }
    
    // return r_global - f_ext_global;
// }


// ///// MATRIX FREE STYLE PRECONDITIONED
// void SolveImplicitStep() {
    // // FORGE-style parameters
    // double beta = 0.25;  // Mass scaling factor
    // double dt = 1.0;     // Pseudo-time step
    
    // Matrix u(m_node_count * m_dim, 1); 
    // Matrix v(m_node_count * m_dim, 1);

    // for (int iter = 0; iter < max_iter; iter++) {
        // // 1. Compute residual
        // Matrix r = ComputeGlobalResidual();
        
        // // 2. FORGE-style mass preconditioning
        // Matrix z(m_node_count * m_dim, 1);
        // for (int n = 0; n < m_node_count; n++) {
            // Matrix r_n = gatherNodalResidual(n, r);
            // z_n = m_mdiag[n].Inv() * r_n / (beta * dt * dt);
            // scatter(z, z_n, n);
        // }
        
        // // 3. Update velocities and positions
        // v += z;
        // u += dt * v;
        
        // // 4. Check convergence
        // if (r.Norm() < tolerance) break;
    // }
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

  double rho_b = 0.818200;  // DEFAULT SPECTRAL RADIUS
  

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

  // if (m_dim == 3 && m_nodxelem ==4){
  // //Replaces PREVIOUS, INSTEAD MASS APPROACH, BUT STILL TO WORK FOR HEXAS
  // cout << "Calc tetra vol"<<endl;
    CalcNodalVol(); //To calc nodal mass
    CalcNodalMassFromVol(); //Repla

  for(int n=0;n<m_node_count;n++)
    for (int d=0;d<m_dim;d++)
      ut_prev[m_dim*n+d]=0.0;
  
  
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
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  while (Time < end_t) {
      
  ////// OR TIME
  if (step_count % 100 == 0){
    printf("Step %d, Time %f, End Time: %.4e, Step Time %.4e\n",step_count, Time, end_t, dt);  
    timer.click();
    //std::cout << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    std::cout << "Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    //std::cout << "CPU Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  }
  
  if (step_count % 10 == 0){
    //cout << "Calc ExtFace Areas"<<endl;
    CalcExtFaceAreas();
    //cout << "Done"<<endl;
    }
  
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

  if (!m_fixed_dt){
    double mat_cs = sqrt(mat[0]->Elastic().BulkMod()/rho[0]);
    calcMinEdgeLength();
    double minl = getMinLength();
    dt = m_cfl_factor*minl/(mat_cs);
  }
  
  //printf("Prediction ----------------\n");
  #if CUDA_BUILD
  N = getNodeCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  #endif
  #if CUDA_BUILD

  #else

  #endif
  
  // !!! PREDICTION PHASE
  // u = dt * (nod%v + (0.5d0 - beta) * dt * prev_a)
  // !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
  // nod%v = nod%v + (1.0d0-gamma)* dt * prev_a
  // nod%a = 0.0d0
  
  // call impose_bcv !!!REINFORCE VELOCITY BC
  // TO NOT INTERFERE WITH DIFF THREADS AND DIMENSIONS
  
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
  //cout <<"Done."<<endl;
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


  calcArtificialViscosity(); //Added to Sigma
  
  
  // Newton-Raphson loop
  double tolerance = 1e-6;
  int max_iter = 10;
  bool converged = false;

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
          
          //double beta = 0.25;
          // // Add mass scaling for stability (FORGE does this)
          // for (int i = 0; i < m_nodxelem*m_dim; i++) {
              // //////    K[i][i] += M_diag[i] / (beta * dt * dt); // beta = 0.25 typically
              // K.Set(i,i,K.getVal(i,i)+ m_mdiag[getElemNode(e,i)] / (beta * dt * dt)); // beta = 0.25 typically
          // }
          

          solver->assembleElement(e, K);
          solver->assembleResidual(e,fint);//SHOULD BE NEGATIVE!
          
          //for (int i=0;i<m_nodxelem;i++)
            

          
          //r_global =r_global + MatMul(K,WHICH MATRIX?);

          cout << "Delta U "<<endl;
          //delta_u_e.Print();
          
      
      } // end par element loop

      m_solver->applyDirichletBCs();
      cout << "Solving system"<<endl;      
      solver->finalizeAssembly();

      m_solver->Solve();
    
    // Update displacements and check convergence
    double max_residual = 0.0;
    for (int n = 0; n < m_node_count; n++) {
        for (int d = 0; d < m_dim; d++) {
            int idx = n * m_dim + d;
            double du = m_solver->getU(n,d);
            
            // Update displacement and position
            //m_solver->addToU(n, d, du);
            u[idx]+=du;x[idx]+=du;
           
            // Track maximum residual
            max_residual = std::max(max_residual, std::abs(du));
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


  

  
  
  #ifdef CUDA_BUILD  
  N = getNodeCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  // UpdateCorrectionPosKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	// cudaDeviceSynchronize();   psfield
  #else

  //IMPLICIT DIFFERENCE
  // Simple position update based on velocities
  for (int i = 0; i < m_node_count * m_dim; i++) {
      x[i] += v[i] * dt;
      
  }

  #endif  

  if (contact){
    double f =1.0;

    if(Time < RAMP_FRACTION*end_t) {
        f = pow(Time/(RAMP_FRACTION*end_t), 0.5);  // Square root for smoother start
    } else {
        f = 1.0;
    }
    for (int n=0;n<trimesh->nodecount;n++){
      trimesh->node_v[n] = f*m_v_orig[n];
      }
      //cout << "Node 0 v"<<(trimesh->node_v[0]).z<<endl;
    }
    
  
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



// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// MAIN SOLVER LOOP (FORGE-STYLE) //////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cout << "Main Loop----" << endl;

// // FORGE-STYLE: Add diagonal mass matrix (lumped)
// Matrix M(m_node_count * m_dim, m_node_count * m_dim);
// M.SetZero();
// for (int e = 0; e < m_elem_count; e++) {
    // for (int n = 0; n < m_nodxelem; n++) {
        // int node = getElemNode(e, n);
        // for (int d = 0; d < m_dim; d++) {
            // int dof = node * m_dim + d;
            // M.Set(dof, dof, M.getVal(dof, dof) + vol[e] * mat[e]->density() / m_nodxelem);
        // }
    // }
// }

// while (Time < end_t) {
    // if (step_count % 10000 == 0)
        // printf("Step %d, Time %f\n", step_count, Time);  

    // // FORGE-STYLE: Store previous configuration
    // for (int n = 0; n < m_node_count; n++)
        // x_old[n] = x[n];  

    // // FORGE-STYLE: Solve for velocities instead of displacements
    // Matrix v(m_node_count * m_dim, 1); // Velocity vector
    // v.SetZero();

    // // Element loop to compute residual
    // for (int e = 0; e < m_elem_count; e++) {
        // // Compute F, strain, stress as before...
        // // ... [keep your existing strain/stress calculation] ...

        // // Compute internal force (same as before)
        // Matrix fint_elem = MatMul(B.getTranspose(), stress_voigt) * vol[e];

        // // Scatter to global residual
        // for (int a = 0; a < m_nodxelem; a++) {
            // int node = getElemNode(e, a);
            // for (int d = 0; d < m_dim; d++) {
                // int dof = node * m_dim + d;
                // v.Set(dof, 0, v.getVal(dof, 0) + fint_elem.getVal(a * m_dim + d, 0));
            // }
        // }
    // }

    // // FORGE-STYLE: Solve M*v = f_ext - f_int (simplified)
    // for (int i = 0; i < m_node_count * m_dim; i++) {
        // if (M.getVal(i, i) > 0) {
            // v.Set(i, 0, (f_ext[i] - v.getVal(i, 0)) / M.getVal(i, i));
        // }
    // }

    // // Apply velocity BCs (replace displacement BCs)
    // for (int dim = 0; dim < m_dim; dim++) {
        // for (int n = 0; n < bc_count[dim]; n++) {
            // int node = (dim == 0) ? bcx_nod[n] : ((dim == 1) ? bcy_nod[n] : bcz_nod[n]);
            // int dof = node * m_dim + dim;
            // v.Set(dof, 0, bc_value / dt); // Convert displacement BC to velocity
        // }
    // }

    // // Update positions
    // for (int i = 0; i < m_node_count * m_dim; i++) {
        // x[i] += v.getVal(i, 0) * dt;
    // }

    // // Optional: Check velocity convergence
    // double max_v = 0.0;
    // for (int i = 0; i < m_node_count * m_dim; i++) {
        // max_v = max(max_v, abs(v.getVal(i, 0)));
    // }
    // if (max_v < velocity_tol) break;

    // step_count++;
    // Time += dt;
// }
