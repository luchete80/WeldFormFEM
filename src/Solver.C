#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

using namespace std;

namespace MetFEM{

	void host_ Domain_d::SolveChungHulbert(){

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

/*  cout << "Initializing bcs" << endl;
    for (int e=0;e<m_elem_count;e++)
      for (int n=0;n<m_nodxelem;n++)
        for (int d=0;d<3;d++)
          v[3*m_elnod[n]+d]=0.0;*/
  
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
  cout << "done"<<endl;
  
  double rho_b = 0.8182;  // DEFAULT SPECTRAL RADIUS
  
  m_alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b);
  m_beta  = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b));
  m_gamma = 1.5 - m_alpha;

  printf ("alpha %f\n", m_alpha);
  printf ("beta %f\n",  m_beta);
  printf ("gamma %f\n", m_gamma);
  
  //cout << "Calculating derivatives..."<<endl;
	#if CUDA_BUILD
	calcElemJAndDerivKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  //cout << "Calculating Volume..."<<endl;
  calcElemInitialVolKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();   
  #else
  calcElemJAndDerivatives();
  CalcElemInitialVol(); //ALSO CALC VOL

  CalcElemVol();
  calcElemDensity();
  calcElemMassMat(); 
  assemblyMassMatrix();  
  
  #endif
	//cout << "Done. "<<endl;


      
  
  Time = 0.0;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Main Loop----"<<endl;
  while (Time < end_t) {
    
  printf("Time %f\n",Time);  
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
 
  //ELEMENT PARALLEL

  #ifdef CUDA_BUILD
  N = getElemCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;  
  
  cout << "calc deriv"<<endl;
	calcElemJAndDerivKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  
  
	calcElemVolKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();   

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
  cout << "dt "<<dt<<endl;
  calcStressStrainKernel<<<blocksPerGrid,threadsPerBlock>>>(this, dt);
  cudaDeviceSynchronize();

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
  //cout << "calc deriv"<<endl;
  calcElemJAndDerivatives();
  //cout << "done "<<endl;
  CalcElemVol();
  //cout << "done "<<endl;
  calcElemStrainRates();
  calcElemDensity();
  calcElemPressure();
  CalcStressStrain(dt);
  calcElemForces();
  calcElemHourglassForces();
  calcElemMassMat(); 
  assemblyMassMatrix();  
  //cout << "Assemblying"<<endl;
  assemblyForces(); 
  //cout << "Done "<<endl;
  calcAccel();
  #endif
  
  ImposeBCAAllDim(); //FOR BOTH GPU AND CPU
  
  N = getNodeCount();
  //printf("Correction\n");	
  #ifdef CUDA_BUILD
  
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  UpdateCorrectionAccVelKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  
  #else
    
  UpdateCorrectionAccVel();
  
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
	cudaDeviceSynchronize();   
  #else
  UpdateCorrectionPos();
  #endif  
  
  // call impose_bcv !!!REINFORCE VELOCITY BC

  // !u = u + beta * nod%v * dt
  // u = u + beta * dt * dt * nod%a   
  // nod%u = nod%u + u
  // nod%x = nod%x + u
  
  // !call AverageData(elem%rho(:,1),nod%rho(:))  
  // prev_a = nod%a
  // time = time + dt
    
  Time += dt;
  
  }// WHILE LOOP


  
  #ifdef CUDA_BUILD
  printf("DISPLACEMENTS\n");
  memcpy__tohost_t(this->u_h,this->u,sizeof(double) * this->m_node_count*3);
  for (int i=0;i<m_node_count;i++)
    printf("%.6e %.6e %.6e\n", this->u_h[i*3+0],this->u_h[i*3+1],this->u_h[i*3+2]);
  printVecKernel<<<1,1 >>>(this, this->u);
	cudaDeviceSynchronize(); 



  memcpy__tohost_t(this->x_h,this->x,sizeof(double) * this->m_node_count*3);

  printf("VELOCITIES\n");
  printVecKernel<<<1,1 >>>(this, this->v);
	cudaDeviceSynchronize(); 
  
  printf("ACCEL\n");
  printVecKernel<<<1,1 >>>(this, this->a);
	cudaDeviceSynchronize();   

  printf("FORCES\n");
  printVecKernel<<<1,1 >>>(this, this->m_fi);
	cudaDeviceSynchronize(); 
  
  printf("STRESSES\n");
  printVecKernel<<<1,1 >>>(this, this->m_sigma);
	cudaDeviceSynchronize();   

  printf("SHEAR STRESSES\n");
  printVecKernel<<<1,1 >>>(this, this->m_tau);
	cudaDeviceSynchronize();
  #else
  calcElemStrainRates();
  printf("DISPLACEMENTS\n");
  printVec(this->u);   

  printf("VELOCITIES\n");
  printVec( this->v);

  printf("ACCEL\n");
	printVec(this->a); 
  
  printf("FORCES\n");
  printVec(this->m_fi);

  printf("STRESSES\n");
  printSymmTens(this->m_sigma);

  printf("SHEAR STRESS\n");
  printSymmTens(this->m_tau);

  printf("STRAIN RATES\n");
  printSymmTens(this->m_str_rate);
  
  printf("ROT RATES\n");
  printSymmTens(this->m_rot_rate);
  
  
  #endif

  VTKWriter writer(this, "out.vtu");
  writer.writeFile();
  
  }//SOLVE
	
};
