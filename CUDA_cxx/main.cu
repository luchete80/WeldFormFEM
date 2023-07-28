#include "Domain.cuh"

#include <iostream>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
//https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void report_gpu_mem()
{
    size_t free, total;
    cudaMemGetInfo(&free, &total);
    std::cout << "Free = " << free << " Total = " << total <<std::endl;
}

using namespace std;

int main() {
  Domain *dom;
  report_gpu_mem();
	gpuErrchk(cudaMallocManaged(&dom, sizeof(Domain)) );
	report_gpu_mem();
  
  double3 V;
  double Lx,Ly,Lz,r;
  Lx = Ly = Lz = 0.1;
  r = 0.05;
  V = make_double3(0.,0.,0.);
  dom->dim = 3;
  
  dom->AddBoxLength(0, V, Lx, Ly, Lz, r);
  
  cout << "Program end."<<endl;
}