#ifndef _CUDAUTILS_H_
#define _CUDAUTILS_H_

#ifdef  CUDA_BUILD
#include "cuda.h"
#include <iostream>

using namespace std;
inline void report_gpu_mem_()
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

  inline  __device__ void double3_Ptr(const double3 v, double *p, const int i){p[3*i]=v.x;p[3*i+1]=v.y;p[3*i+2]=v.z;}
  inline  __device__ double3 Ptr_double3(double *p, const int i){return make_double3(p[3*i],p[3*i+1],p[3*i+2]);}
  inline  __device__ double2 Ptr_double2(double *p, const int i){return make_double2(p[3*i],p[3*i+1]);}
#else
  #include "double3_c.h"
  //inline void Vec3D_Ptr(const Vec3D v, double *p, const int i){p[3*i]=v(0);p[3*i+1]=v(1);p[3*i+2]=v(2);}
  //inline Vec3D Ptr_Vec3D(double *p, const int i){return Vec3D(p[3*i],p[3*i+1],p[3*i+2]);}
  inline void double3_Ptr(const double3 v, double *p, const int i){p[3*i]=v.x;p[3*i+1]=v.y;p[3*i+2]=v.z;}
  inline double3 Ptr_double3(double *p, const int i){return make_double3(p[3*i],p[3*i+1],p[3*i+2]);}
  
  //// TEMPLATIZE THIS
  inline double2 Ptr_double2(double *p, const int i){return make_double2(p[2*i],p[2*i+1]);}
  inline void double2_Ptr(const double2 v, double *p, const int i){p[2*i]=v.x;p[2*i+1]=v.y;}
#endif

#endif
