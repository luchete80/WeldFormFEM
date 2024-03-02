#include "Domain_d.h"
#include <iostream>

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
  AssignMatAddressKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
  cudaDeviceSynchronize();
  #else 
  AssignMatAddress();
  #endif

  for (int d=0;d<m_dim;d++){
    N = bc_count[d];
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    
    #ifdef CUDA_BUILD
    ImposeBCVKernel<<<blocksPerGrid,threadsPerBlock >>>(this, d);
    cudaDeviceSynchronize();
    #else
      ImposeBCV(d);
    #endif
  }
  
  
  double rho_b = 0.8182;  // DEFAULT SPECTRAL RADIUS
  
  m_alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b);
  m_beta  = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b));
  m_gamma = 1.5 - m_alpha;

  printf ("alpha %f\n", m_alpha);
  printf ("beta %f\n",  m_beta);
  printf ("gamma %f\n", m_gamma);
  
  
  cout << "Calculating derivatives..."<<endl;
	#if CUDA_BUILD
	calcElemJAndDerivKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  cout << "Calculating Volume..."<<endl;
  calcElemInitialVolKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();   
  #else
  calcElemJAndDerivatives();
  CalcElemInitialVol();
  #endif
	cout << "Done. "<<endl;
  
  Time = 0.0;
  
  while (Time < end_t) {

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
    N = bc_count[d];
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    
    #ifdef CUDA_BUILD
    ImposeBCVKernel<<<blocksPerGrid,threadsPerBlock >>>(this, d);
    cudaDeviceSynchronize();
    #else
      ImposeBCV(d);
    #endif
  }
  cout <<"Done."<<endl;

  N = getElemCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;  
  
  //ELEMENT PARALLEL

  #ifdef CUDA_BUILD
	calcElemJAndDerivKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 

	calcElemVolKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize();   

  calcElemStrainsKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 

  calcElemDensityKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 

  calcElemPressureKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
  cudaDeviceSynchronize();   
  #else
  calcElemJAndDerivatives();
  #endif
  // calcElemForcesKernel<<<blocksPerGrid,threadsPerBlock>>>(this);
  // cudaDeviceSynchronize();   



  N = getNodeCount();
  blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

  #ifdef CUDA_BUILD
  UpdateCorrectionKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); 
  
  #else
    
  UpdateCorrection();
  
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
  
  }//SOLVE
	
};