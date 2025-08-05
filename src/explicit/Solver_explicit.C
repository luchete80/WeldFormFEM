#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

#include "Mesh.h"
#include "WallTimer.h"

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>  // ofstream
#include <iomanip>

#ifdef BUILD_REMESH
#include "ReMesher.h"
#endif

using namespace std;



std::ostringstream m_oss;
std::string m_fname;


namespace MetFEM{

void host_ Domain_d::SolveChungHulbert(){
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
  
  m_alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b);
  m_beta  = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b));
  m_gamma = 1.5 - m_alpha;

  printf ("alpha %.10e\n", m_alpha);
  printf ("beta %.10e\n",  m_beta);
  printf ("gamma %.10e\n", m_gamma);
  
  ostringstream oss_out;
  
  if (m_press_algorithm==0){
    oss_out << "alpha_free "<<m_stab.alpha_free <<endl;
    oss_out << "alpha_contact "<<m_stab.alpha_contact <<endl;
    oss_out << "hg_coeff_free "<<m_stab.hg_coeff_free <<endl;
    oss_out << "hg_coeff_contact "<<m_stab.hg_coeff_contact <<endl;
    oss_out << "av_coeff_div "<<m_stab.av_coeff_div <<endl;
    oss_out << "av_coeff_bulk "<<m_stab.av_coeff_bulk <<endl;
    oss_out << "log_factor "<<m_stab.log_factor <<endl;
    oss_out << "pspg_scale "<<m_stab.pspg_scale <<endl;
    oss_out << "p_pspg_bulkfac "<<m_stab.p_pspg_bulkfac <<endl;
    
    cout << oss_out.str();
    out_file << oss_out.str();
    out_file.flush();
  }
  
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

  if (m_dim == 2 && m_domtype == _Axi_Symm_)
    Calc_Element_Radius();
  
  CalcElemInitialVol(); //ALSO CALC VOL
  
  CalcElemVol();
  printf("calc dens\n");
  calcElemDensity();

  // if (m_dim == 3 && m_nodxelem ==4){
  // //Replaces PREVIOUS, INSTEAD MASS APPROACH, BUT STILL TO WORK FOR HEXAS
  // cout << "Calc tetra vol"<<endl;
    CalcNodalVol(); //To calc nodal mass
    CalcNodalMassFromVol(); //Repla
  //ONLY FOR 4 elem cube
  double totmass = 0.0;
  double totvol= 0.0;
  for(int n=0;n<m_node_count;n++)
   totmass +=m_mdiag[n];
 
   for(int e=0;e<m_elem_count;e++)
   totvol +=vol[e];
    
  
  cout << "TOTAL MASS: "<<totmass<<endl;
  cout << "TOTAL VOL: "<<totvol<<endl;  

    // calcElemMassMat();
    // assemblyMassMatrix();    
    
    // }
  //CONTACT:
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

  Time = 0.0;
  int step_count = 0;
  double tout = 0;
  
  bool remesh_ = false;
  int remesh_count = 0;
  const double RAMP_FRACTION = 1.0e-2;  // 0.1% of total time instead of 1%
  of << "t,f,area"<<endl;
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
  if (remesh_){
      cout << "New dt: "<< dt<<endl;
  }
  
  //printf("Prediction ----------------\n");
  #if CUDA_BUILD
  N = getNodeCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  #endif
  #if CUDA_BUILD
  UpdatePredictionKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  #else
  UpdatePrediction();  
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
    if (m_dim == 2 && m_domtype == _Axi_Symm_)
      Calc_Element_Radius();
    
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
  
  calcAccelKernel<<<blocksPerGrid,threadsPerBlock>>>(this);
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
  //smoothPressureField(0.2);
  //calcElemPressure_Hybrid_VolHG();
  //calcElemPressureANP_Nodal_HG();
  
  //calcElemPressureFromJ();
  CalcStressStrain(dt);
  calcArtificialViscosity(); //Added to Sigma
  
  calcElemForces();
  calcElemHourglassForces();
  
  if (contact)
//    CalcContactForcesWang();
    CalcContactForces();

  //if (m_dim == 3 && m_nodxelem ==4){
  //Replaces PREVIOUS, INSTEAD MASS APPROACH, BUT STILL TO WORK FOR HEXAS
    //THIS IS NOT WORKING
    //CalcNodalVol(); //To calc nodal mass
    //CalcNodalMassFromVol(); //Repla
  //} else{

    //calcElemMassMat();
    //assemblyMassMatrix();    
    
    //}
     
  //calcElemMassMat(); 
  //assemblyMassMatrix();  
  
  assemblyForces(); 
  //ApplyGlobalSprings();

  calcAccel();
  
  #endif
  
  ImposeBCAAllDim(); //FOR BOTH GPU AND CPU
  
  N = getNodeCount();
  //printf("Correction\n");	
  #ifdef CUDA_BUILD
  
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    //printf("Upate Correction\n");
  UpdateCorrectionAccVelKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  
  #else
  UpdateCorrectionAccVel();
  
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

	// !$omp parallel do num_threads(Nproc) private (n)  
	// do n=1,node_count
		// ! do d=1,dim
			// ! nod%a(n,d) =  (fext_glob(n,d)-rint_glob(n,d))/mdiag(n) 
		// ! end do 
		// nod%a(n,:) =  (fext_glob(n,:)-rint_glob(n,:))/mdiag(n) 
	// end do
	// !$omp end parallel do
	
  // call impose_bca

  
	// !$omp parallel do num_threads(Nproc) private (n)
  // do n=1,node_count
		// nod%a(n,:) = nod%a(n,:) - alpha * prev_a(n,:)
		// nod%a(n,:) = nod%a(n,:) / (1.0d0 - alpha)
		// nod%v(n,:) = nod%v(n,:) + gamma * dt * nod%a (n,:)  
	// end do
	// !$omp end parallel do
  

  ImposeBCVAllDim();
  
  #ifdef CUDA_BUILD  
  N = getNodeCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  UpdateCorrectionPosKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();   psfield
  #else
  UpdateCorrectionPos();
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

  // call impose_bcv !!!REINFORCE VELOCITY BC

  // !u = u + beta * nod%v * dt
  // u = u + beta * dt * dt * nod%a   
  // nod%u = nod%u + u
  // nod%x = nod%x + u
  
  // !call AverageData(elem%rho(:,1),nod%rho(:))  
  // prev_a = nod%a
  // time = time + dt
 
  if (Time>=tout){
    string outfname = "out_" + std::to_string(Time) + ".vtk";
    timer.click();

    ostringstream oss_out;
    oss_out << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    oss_out << "CPU Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    oss_out << "Plastic Strain energy "<<m_pl_energy<<endl;
    //printf("Reaction Forces\n");
    of <<std::scientific<<std::setprecision(6)<< Time ;
    //for (int m=0;m<trimesh->mesh_count;m++){
      //printf("Surf Id %d %.4e\n",m, norm(trimesh->react_force[m]));
    //printf("Surf Id %d %.4e\n",m, trimesh->react_p_force[m]);
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
  //ReMesher remesh(this);
  
  //remesh.Generate_mmg();
  //remesh.m_type = MMG;
  #endif
  //////////////////////////////////////
  

  #ifdef CUDA_BUILD
  cudaMemcpy(x_h, x, 3*sizeof(double) * m_node_count, cudaMemcpyDeviceToHost);		
  
  //printf("X %.6e\n", x_h[0]); //CRASHES

/*
  printf("DISPLACEMENTS\n");
  printVecKernel<<<1,1 >>>(this, this->u);
	cudaDeviceSynchronize(); 

  printf("VELOCITIES\n");
  printVecKernel<<<1,1 >>>(this, this->v);
	cudaDeviceSynchronize(); 

  printf("ACCEL\n");
  printVecKernel<<<1,1 >>>(this, this->a);
	cudaDeviceSynchronize(); 

  printf("FORCES\n");
  printVecKernel<<<1,1 >>>(this, this->m_fi);
	cudaDeviceSynchronize(); 
*/

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
