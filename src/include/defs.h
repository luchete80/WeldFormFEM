#ifndef _DEFS_H_
#define _DEFS_H_

//TODO: TEST IT AND COMPARE WITH VEC3D...
#include "Vec3D.h"



#ifdef  CUDA_BUILD
#include "utils.h"

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
#define memcpy__tohost_t(dest,src,size)   cudaMemcpy(dest,src,size,cudaMemcpyDeviceToHost)

#define par_loop(n,upto)   blocksPerGrid = (upto + threadsPerBlock - 1) / threadsPerBlock; int n=threadIdx.x+blockDim.x*blockIdx.x;if(n<upto)

//--------------------------------------------------------------
#else


#include "double3_c.h"
#include "utils.h" //AFTER DEF

#include "Vec3D.h"
#define __global__ 
#define spec_

#define vector_t double3
#define dev_t
#define make_vector_t make_double3

#include <cmath> //POW
#define host_ 

#define vector_t_Ptr  double3_Ptr
#define Ptr_vector_t  Ptr_double3
#define Ptr_vector2   Ptr_double2
//#define V_.x   V_[0]
// inline void Vec3D_Ptr(const Vec3D v, double *p, const int i);
// inline Vec3D Ptr_Vec3D(double *p, const int i);
// #define Ptr_vector_t  Ptr_Vec3D
// #define vector_t_Ptr  Vec3D_Ptr

#define malloc_t(x,t,y)           x=(t*)malloc(y * sizeof(t))
#define memcpy_t(dest,src,size)   memcpy(dest,src,size)     
#define OMP_PARA_INTERNAL _Pragma("omp parallel for schedule (static) num_threads(Nproc)")

#define par_loop(n,upto)   OMP_PARA_INTERNAL\
                            for(int n=0;n<upto;n++)

//#define par_loop(n,upto)  for(int n=0;n<upto;n++)
#endif


// FOR CUDA AND CPU
//#define make_vector_t Vec3D
//#define make_vector_t make_double3
//#if 
//inline make_vector_t(){}

#endif


