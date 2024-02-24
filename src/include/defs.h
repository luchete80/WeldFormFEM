#ifndef _DEFS_H_
#define _DEFS_H_

#ifdef  CUDA_BUILD

#define vector_t double3
#define spec_ __host__ __device__
#define dev_t __device__
#define host_ __host__
#else

#include "Vec3D.h"

#define spec_
#define vector_t Vec3D
#define dev_t
#include <cmath> //POW
#define host_ 
#endif

#endif
