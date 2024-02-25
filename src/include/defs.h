#ifndef _DEFS_H_
#define _DEFS_H_

#include "Vec3D.h"
#include "utils.h"

#ifdef  CUDA_BUILD

#define vector_t double3
#define spec_ __host__ __device__
#define dev_t __device__
#define host_ __host__
#define make_vector_t make_double3
#define vector_t_Ptr  double3_Ptr
#define Ptr_vector_t  Ptr_double3
#define loop (n,upto)   n=threadIdx.x+blockDim.x*blockIdx.x;if(n<upto)

#define malloc_t(x,t,y)           cudaMalloc((void **)&x, y * sizeof (t))
#define memcpy_t(dest,src,size)   cudaMemcpy(dest,src,size,cudaMemcpyHostToDevice)

#define par_loop(n,upto)    int n=threadIdx.x+blockDim.x*blockIdx.x;if(n<upto)

#else

#include "Vec3D.h"
#define __global__ 
#define spec_
#define vector_t Vec3D
#define dev_t
#include <cmath> //POW
#define host_ 
#define make_vector_t Vec3D
inline void Vec3D_Ptr(const Vec3D v, double *p, const int i);
#define vector_t_Ptr  Vec3D_Ptr

#define malloc_t(x,t,y)           x=(t*)malloc(y * sizeof(t))
#define memcpy_t(dest,src,size)   memcpy(dest,src,size)     
#define par_loop(n,upto)   int n;for(int n=0;n<upto;n++)
#endif

#endif
