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
  
  
  //////////// TEST ; TO DELETE /////
  bool elem_test[m_elem_count];
  for (int e=0;e<m_elem_count;e++)elem_test[e]=false;


   //double pxa_el[m_elem_count];
  double zmax = 0.0;
      for (int i=0;i<m_node_count;i++)
        if (getNodePos3(i).z>zmax)
          zmax = getNodePos3(i).z;
        
    int ecount = 0;
   double area = 0.0;
   

  double cfsum = 0.0;
  for (int i=0;i<m_node_count;i++){
    if (abs(getNodePos3(i).z-zmax)<1.0e-5){
      for (int ne=0; ne<m_nodel_count[i];ne++) {
        int e   = m_nodel     [m_nodel_offset[i]+ne]; //Element
        if (!elem_test[e]){
          //pxa_el[e]+=p[e]*m_elem_area[e];
          elem_test[e]=true;
          cfsum += p[e]*m_elem_area[e];
          area+=m_elem_area[e];
          //cout << "element "<<e<<endl;
          ecount++;
        }
      }//nodel
    }//if inside
  }//node
  cout << "UPPER ELEMENT COUNT"<<ecount<<endl;

  ///////////////// END TEST ///
  
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
    oss_out << "log_factor     "<<m_stab.log_factor <<endl;
    oss_out << "pspg_scale     "<<m_stab.pspg_scale <<endl;
    oss_out << "p_pspg_bulkfac "<<m_stab.p_pspg_bulkfac <<endl;
    oss_out << "J_min          "<<m_stab.J_min <<endl;
    oss_out << "hg_forces      "<<m_stab.hg_forces <<endl;
            
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
  of << "t,f,fc,area,Eint,Ekin"<<endl;
  int last_step_remesh=-1;
  double Ekin, Eint;
  double Ekin_old;
  Ekin = Eint = 0.0;
  double dEkin,dEint;
  bool decrease_dt = false;
  double max_vel = 0.0;
  int wup_step_count = 0;
  bool transition = false;
  int trans_step_count = 0;
  int end_wup_step;
  double r_damp;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double max_vprev=0.0;
  double s_wup = 1.0; //WARM UP STEP
  while (Time < end_t) {
      
  ////// OR TIME
  if (step_count % 100 == 0){
    printf("Step %d, Time %f, End Time: %.4e, Step Time %.4e\n",step_count, Time, end_t, dt);  
    cout << "Ekin: "<<Ekin <<", Eint "<<Eint<<endl; 
    cout << "Max V"<<max_vel<<endl;
    timer.click();
    //std::cout << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    std::cout << "Overall Time" << timer.elapsedSinceStart() << " seconds\n";
    //std::cout << "CPU Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  }
  

  
  //~ if (step_count % 50 == 0)
    //~ SearchExtNodes(); //TODO: CALCULATE ONLY AREA, NOT SEARCH AGAIN AREAS

  double initial_dt;
  /////AFTER J AND DERIVATIVES
  if ( (step_count - last_step_remesh) % m_remesh_interval == 0 && step_count  >0 && remesh_count < m_remesh_max_count)
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
      
      //smoothFieldLaplacian(v,3);
      s = "out_remesh_smooth"+std::to_string(step_count)+".vtk";
      VTKWriter writer4(this, s.c_str());
      writer4.writeFile();
      Ekin_old = Ekin;
      #endif
      remesh_count++;
      last_step_remesh = step_count;
      s_wup= 0.0;
      max_vprev = 0.0;
      wup_step_count = 0;
      transition = false;
      r_damp = 1.0;
      }
      //#########################################################
  //////////////////////////// IF REMESH

  }

  if (m_dim > 2){
    if (step_count % 10 == 0 || remesh_){
    //cout << "Calc ExtFace Areas"<<endl;
    CalcExtFaceAreas();
    //cout << "Done"<<endl;
    }
  }

  max_vel=0.0;
  if (!m_fixed_dt){
    double mat_cs = sqrt(mat[0]->Elastic().BulkMod()/rho[0]);
    calcMinEdgeLength();
    double minl = getMinLength();
    for (int i=0;i<m_node_count;i++){
      vector_t v = getVelVec(i);
      if(norm(v)>max_vel)
        max_vel = norm(v);
    
    double prev_dt = dt;
    dt = m_cfl_factor*minl/(mat_cs+max_vel);
    if (remesh_)
      initial_dt = dt;
    if (dt>10.0*prev_dt) cout << "ERROR: DT "<<dt<<endl;
      
    }

    // if (dt_new > 1.0e-20)
      // dt = dt_new;
    // else{
      // ////DIVERGENCE
      // cout << "ERROR DT ZERO!, GETTING PREVIOUS---"<<endl;
      
    // }
  }
  
  

  #ifdef BUILD_REMESH  
  if (s_wup < 1.0){
    //s_wup = double(step_count-last_step_remesh)/double(m_filter_params.warmup_steps);

    //LINEAR
    //~ if (!decrease_dt) //maintain dt instead, to ensure it will converge
      //~ s_wup += 1.0/double(m_filter_params.warmup_steps);

    //CUBIC EASE OUT
    if (!decrease_dt){
      double s_norm = double(step_count - last_step_remesh) / double(m_filter_params.warmup_steps);
      s_norm = std::min(1.0, s_norm);
      s_wup = 1.0 - pow(1.0 - s_norm, 3); // cubic ease-out
    }
    
    if (wup_step_count % 10 == 0){
      
      string s = "out_wup_"+std::to_string(wup_step_count) + std::to_string(step_count)+".vtk";
      VTKWriter writer(this, s.c_str());
      writer.writeFile();
    }
      cout << "s warmup: "<<s_wup<<endl;
    // 1. Calcular dt_base independientemente del dt anterior
    double dt_base = dt * 0.8 * pow(s_wup+(1.0/double(m_filter_params.warmup_steps)), 2.0);
    
    // 2. Aplicar límites físicos
    double dt_CFL = dt;
    dt = std::min(dt_base, dt_CFL);
    
    // 3. Nunca permitir dt=0
    dt = std::max(dt, 1e-12);
    
    cout << "s warmup: " << s_wup << endl;
    cout << "New dt: " << dt << endl;
    cout << "Max vel: " << max_vel << endl;

    cout << "dEkin: "<<dEkin <<", dEint: "<<dEint<<endl;
    wup_step_count++;
    
      // std::string s = "out_wup_"+std::to_string(s_wup)+".vtk";
      // VTKWriter writer3(this, s.c_str());
      // writer3.writeFile();
  } else {
    
        if (!transition && wup_step_count > 0) {
          /*
          if (max_vel>10.0){  //10 MESH SPEED
          
            s_wup = r_damp;
            cout << "CORRECTED WARMUP according to damping"<<endl;
          }else {
        */
          transition = true;
          cout << "Transition phase. Warmup completed in " << wup_step_count<<" steps."<<endl;
          trans_step_count = 0;
          end_wup_step = step_count;
          //smoothFieldLaplacian(v,3);
          string s = "out_remesh_wup_bef_trans"+std::to_string(step_count)+".vtk";
          VTKWriter writer(this, s.c_str());
          writer.writeFile();
      //}
        
    }
    
    }
  
  // if (decrease_dt){
    // //m_filter_params.warmup_steps *=2;
    // //s_wup = 0.5*1.0/(double(m_filter_params.warmup_steps));
    // s_wup -= 2.0*1.0/double(m_filter_params.warmup_steps);
  // }
  
  
  #endif
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
  
  #ifdef BUILD_REMESH    
  if (remesh_){
  //Maintain mapped v if were mapped with momentum conservation!
    memcpy_t(v,  m_vprev, sizeof(double) * m_node_count * 3); 
  }
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
  //smoothDevStrainRates(0.5);
  #ifdef BUILD_REMESH    
  if (s_wup < 1.0 ){  
    BlendField(s_wup,m_elem_count,6,m_str_rate_prev,m_str_rate);
  }
  #endif
  

  // if (m_dim == 3 && m_nodxelem ==4){
  //calcElemPressureANP_Nodal();
  //calcElemPressureANP();
  // }else
  //if (!remesh_){
  if      (m_press_algorithm == 0)
    calcElemPressure();
  else if (m_press_algorithm == 1)
    calcElemPressureANP();
  //}

  //if (step_count%5 == 0)
    //smoothPressureLaplacian();
  
  calcNodalPressureFromElemental();
  //smoothPressureField(0.2);
  //calcElemPressure_Hybrid_VolHG();
  //calcElemPressureANP_Nodal_HG();
  
  //calcElemPressureFromJ();
  
  CalcStressStrain(dt);
  calcArtificialViscosity(); //Added to Sigma

  #ifdef BUILD_REMESH    
  if (s_wup < 1.0 ){
    BlendStresses(s_wup, 1.5);
    SmoothDeviatoricStress(0.8);
  }
  #endif
  
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
  
  #ifdef BUILD_REMESH    
  // DIVERGENCE DETECTION 
  // COULD BE ALSO WITH KINETIC ENERGY  
  int nc=0;
  bool large_acc = false;
  double maxv = 0.0;
  for (int i=0;i<m_node_count;i++){
      vector_t acc = getAccVec(i);
      vector_t vel = getVelVec(i);
      if(norm(acc)>1.0e10 ){
        //nc++;
        //large_acc = true;
        for (int d=0;d<m_dim;d++){
        
            //postRemeshGlobFilter();
            a[m_dim*i+d] *= 0.01;
        }
      }

      if(norm(vel)>20.0*1.2){
        if (norm(vel)>maxv)
          maxv =norm(vel);
        nc++;

      }

  }
  
  if (nc>0.05*m_node_count)
        large_acc = true;
      
  if (large_acc){ 
    cout << "ERROR, "<< nc <<" nodes with veloc too large "<<", max vel: "<<maxv<<endl;  
    decrease_dt = true;
  } else {
    decrease_dt = false;
    }
    
    
  #endif 
  
  #endif //CUDA_BUILD
  
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
  
  #ifdef BUILD_REMESH    
  if (remesh_){
    //Maintain mapped v if were mapped with momentum conservation!
    memcpy_t(v,  m_vprev, sizeof(double) * m_node_count * 3); 
  }
  
  double v_max = 0.0;

  
  if(s_wup<1.0|| transition){
     r_damp = sqrt( Ekin_old / (Ekin + 1e-30) );    
      for (int n=0;n<m_node_count;n++){ 
      vector_t vel = getVelVec(n);
      if(norm(vel)>v_max ){
        v_max = norm(vel);
      }
      
        for (int d=0;d<m_dim;d++)
          if (r_damp<1.0)
            v[m_dim*n+d] *= r_damp;   // nunca subir v; solo bajar si se disparó    

    CorrectLocalVelocityPeaks();
      }
  }
  

  
  if (contact){
    //if (Time > RAMP_FRACTION*end_t)
    //ApplyGlobalDamping(0.02);
  }
  
  int DAMPING_TRANSITION_STEPS = 100;
  
  if (s_wup < 1.0 || transition){
    if (max_vprev>0.0){ //NOT FIRST WARM UP STEP
      if (v_max>max_vprev*1.1){
        decrease_dt = true;
        cout << "V>1,1PREV_V"<<endl;
      }
    }
    max_vprev = v_max;
    cout << "Max vel before correct peaks and affect with Ekin"<<v_max<<endl;
    cout << "Damping Energy factor: "<<r_damp<<endl;
    //if (Time > RAMP_FRACTION*end_t)
    
    //ApplyGlobalDamping((1.0-s_wup));
    
    //smoothFieldLaplacian(v,3);
    const double ka = 0.5;
    for (int i=0;i<m_node_count;i++)
      for (int d=0;d<m_dim;d++){
        //if(abs(a[m_dim*i+d])>1.0e6)
      
          postRemeshGlobFilter();//Filter again vel: TOO SLOW!!
          a[m_dim*i+d] *= ka*double(step_count-last_step_remesh)/double(m_filter_params.warmup_steps);
          //v[m_dim*i+d] *= (1.0e-2)*double(step_count-last_step_remesh)/double(m_filter_params.warmup_steps);
      }

  }
  #endif //REMESH
  
  if (transition) {
    ///if (step_count - end_wup_step < m_filter_params.trans_step_count){
    if (trans_step_count < m_filter_params.trans_step_count){
        //DO NOT BLEND
        for (int n=0;n<m_node_count;n++){ 
        vector_t vel = getVelVec(n);
        if(norm(vel)>v_max )
          v_max = norm(vel);
        }    
        cout << "Transition step "<< trans_step_count << " / " <<  m_filter_params.trans_step_count <<", vmax "<<v_max<<endl;
        //postRemeshGlobFilter();
        for (int i=0;i<m_node_count;i++)
          for (int d=0;d<m_dim;d++)
            a[m_dim*i+d] *= 0.75;
        trans_step_count++;
        
    } else {
      
      transition = false;
      wup_step_count = false; //To not reactivate 
      cout << "End transition "<<endl;

      //smoothFieldLaplacian(v,3);
      string s = "out_remesh_wup"+std::to_string(step_count)+".vtk";
      VTKWriter writer(this, s.c_str());
      writer.writeFile();

    }
  }
  //ApplyGlobalDamping(0.1);
  #endif 

  ImposeBCVAllDim();
  //DAMPING
  int nc_axis=0;
  if (m_domtype == _Axi_Symm_){
    double xmin = 1000.0;
    for (int i=0;i<getNodeCount();i++)
      if (getPosVec2(i).x < xmin)
        xmin = getPosVec2(i).x;
        
    for (int i=0;i<getNodeCount();i++){
      if (getPosVec2(i).x <= xmin+1.e-6){
        a[m_dim*i] = 0.0;
        v[m_dim*i] = 0.0;
        nc_axis++;
      }
    }//nodes
    //cout << "nc  axis  "<< nc_axis<<endl;
  }

  
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
  double dEvisc;
  computeEnergies(dt,dEkin,dEint,dEvisc);
  Eint+=dEint; Ekin +=dEkin;

 
  if (Time>=tout){
    string outfname = "out_" + std::to_string(Time) + ".vtk";
    timer.click();

    ostringstream oss_out;
    oss_out << "Step Time" << timer.elapsedSinceLastClick() << " seconds\n";
    cout << "step time "<<endl;
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
                                                    trimesh->react_force[0].z<<","<<
                                                    trimesh->cont_area;
    } else{


    /////////////// TEST NO CONTACT
    double cfs=0.0;
    double cfa = 0.0;
    for (int e=0;e<m_elem_count;e++){
      double f = 1.0;
      if (m_nodxelem == 8) f = 2.0;
      if (elem_test[e]){
        cfa += f*m_elem_area[e];
        cfs += m_sigma[6*e+2]*f*m_elem_area[e];
      }
    }
    
    cout << "cfs "<<cfs<<", cfa "<<cfa;
    //// TEST
    of << ", "<<cfs<<","<<cfa;
    
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
      //cout << "Cont Elements"<<ecount<<endl;
      //of <<std::scientific<<std::setprecision(6)<<", "<<cfsum;
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
        out_file << oss_out.str();
        out_file.flush();
    }
    
    
    //}
    of <<","<<Eint<<","<<Ekin;
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
