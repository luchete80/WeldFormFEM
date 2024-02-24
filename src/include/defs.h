#ifndef _DEFS_H_
#define _DEFS_H_

#ifdef  CUDA_BUILD

#define vector_t double3
#define spec_ __host__ __device__
#define dev_t __device__
#define host_ __host__
#define make_vector_t make_double3
#define vector_t_Ptr  double3_Ptr
#define Ptr_vector_t  Ptr_double3

// typename <T>
// void malloc_t (T,size_t size){
	// cudaMalloc((void **)&x, node_count * sizeof (double) * 3);
// }
#define malloc_t(x,t,y) cudaMalloc((void **)&x, y * sizeof (t))
#else

#include "Vec3D.h"

#define spec_
#define vector_t Vec3D
#define dev_t
#include <cmath> //POW
#define host_ 
#define make_vector_t Vec3D

#define malloc_t(x,t,y) x = (t*) malloc(y * sizeof(t))

#endif

#endif
