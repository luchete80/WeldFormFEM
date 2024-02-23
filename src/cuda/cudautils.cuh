#ifndef _CUDAUTILS_H_
#define _CUDAUTILS_H_

#include "cuda.h"

void    report_gpu_mem_();
inline  __device__ void double3_Ptr(const double3 v, double *p, const int i){p[3*i]=v.x;p[3*i+1]=v.x;p[3*i+2]=v.z;}
inline  __device__ double3 Ptr_double3(double *p, const int i){return make_double3(p[3*i],p[3*i+1],p[3*i+2]);}

#endif