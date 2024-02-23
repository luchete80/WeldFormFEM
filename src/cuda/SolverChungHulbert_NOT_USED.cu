#include "Solver.cuh"
#include <iostream>

using namespace std;

namespace MetFEM{


SolverChungHulbert::SolverChungHulbert(Domain_d *d){
	m_dom = d;
}

__host__ void SolverChungHulbert::Solve(){
	
	int N = m_dom->getElemCount();
	m_dom->threadsPerBlock = 256; //Or BlockSize
	//m_dom->threadsPerBlock = 1; //Or BlockSize
	m_dom->blocksPerGrid =				// Or gridsize
	(N + m_dom->threadsPerBlock - 1) / m_dom->threadsPerBlock;
	cout << "Blocks per grid"<<m_dom->blocksPerGrid<<", Threads per block"<< m_dom->threadsPerBlock<<endl;
	
	calcElemJAndDerivKernel<<<m_dom->blocksPerGrid,m_dom->threadsPerBlock >>>(m_dom);
	cudaDeviceSynchronize(); 

  // !!! PREDICTION PHASE
  // u = dt * (nod%v + (0.5d0 - beta) * dt * prev_a)
  // !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
  // nod%v = nod%v + (1.0d0-gamma)* dt * prev_a
  // nod%a = 0.0d0
  
  // call impose_bcv !!!REINFORCE VELOCITY BC
  // TO NOT INTERFERE WITH DIFF THREADS AND DIMENSIONS
  for (int d=0;d<m_dim;d++){
    ImposeBCVKernel<<<m_dom->blocksPerGrid,m_dom->threadsPerBlock >>>(m_dom, d);
    cudaDeviceSynchronize();
  }
  
  // calcElemStrains<<<m_dom->blocksPerGrid,m_dom->threadsPerBlock >>>(m_dom);
	// cudaDeviceSynchronize(); 
    
  // calcElemPressureKernel<<<m_dom->blocksPerGrid,m_dom->threadsPerBlock >>>(m_dom);
  // cudaDeviceSynchronize();   

  // calcElemForcesKernel<<<m_dom->blocksPerGrid,m_dom->threadsPerBlock >>>(m_dom);
  // cudaDeviceSynchronize();   
  
}

}; //Namespace