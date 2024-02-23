#ifndef _DEFS_H_
#define _DEFS_H_

#ifdef  CUDA_BUILD
#define spec_ __host__ __device__
#define dev_ __device__
#else
#define spec_
#define dev_
#include <cmath> //POW
#endif

#endif
