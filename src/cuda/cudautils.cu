#include <iostream>

#include "cudautils.cuh"



using namespace std;
void report_gpu_mem_()
{
  size_t free, total;
  cudaMemGetInfo(&free, &total);
  std::cout << "Free = " << free << " Total = " << total <<std::endl;
  
  float free_m,total_m,used_m;
  

  free_m =(unsigned int)free/1048576.0 ;

  total_m=(unsigned int)total/1048576.0;

  used_m=total_m-free_m;

  //printf ( "  mem free %d .... %f MB mem total %d....%f MB mem used %f MB\n",free_t,free_m,total_t,total_m,used_m);
  
  std::cout << "  mem free "<< free_m <<" MB mem total" << total_m <<" MB mem used "<< used_m<<"MB"<<endl;
    
}