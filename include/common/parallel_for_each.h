/*************************************************************************/
/*  parallel_for_each.h                                          */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                              */
/*  https://www.opensourcemech.com                                                                */
/*                                                                       */
/*  Copyright (c) 2025-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/



#ifndef PARALLEL_FOR_EACH_H
#define PARALLEL_FOR_EACH_H

#include <iterator>
#include <utility>

// Preprocessor flags for enabling CUDA or OpenMP
#if defined(CUDA_BUILD)
#include <cuda_runtime.h>
#endif

#if defined(USE_OPENMP)
#include <omp.h>
#endif

#ifndef __CUDACC__  // If not compiling with NVCC
    #define __host__
    #define __device__
#endif


namespace parallel {

#if defined(CUDA_BUILD)


// CUDA kernel for parallel for_each
//template <class F, class ForwardIt>
//__global__
//void cuda_for_each(F f, ForwardIt first, ForwardIt last) {
//  using difference_type = typename std::iterator_traits<ForwardIt>::difference_type;
//  auto const i = static_cast<difference_type>(threadIdx.x + blockIdx.x * blockDim.x);
 // ForwardIt const it = first + i;
//  if (it < last) f(*it);
//}


// CUDA kernel for parallel for_each specialized for double*
template <class F>
__global__
void cuda_for_each(F f, double* first, double* last) {
    auto const i = threadIdx.x + blockIdx.x * blockDim.x;
    double* it = first + i;
    if (it < last) f(*it);
}


// Helper function for ceiling division
template <class T>
__host__ __device__ inline constexpr
T ceildiv(T a, T b) {
  return (a / b) + ((a % b) ? 1 : 0);
}

#endif

// Unified for_each function
template <class ForwardIt, class UnaryFunction>
void for_each(ForwardIt first, ForwardIt last, UnaryFunction f) {
  auto const n = last - first;
  printf("N THREADS %d \n",n);
  if (n <= 0) return;

#if defined(CUDA_BUILD)
  // GPU execution path
  printf("Executing CUDA\n");
  dim3 const cuda_block(32, 1, 1);
  //printf("grid %s\n"parallel::ceildiv(static_cast<unsigned>(n)) );
  dim3 const cuda_grid(parallel::ceildiv(static_cast<unsigned>(n), cuda_block.x), 1, 1);
  std::cout << "Launching kernel with " << cuda_grid.x << " blocks and " << cuda_block.x << " threads per block." << std::endl;
  std::size_t const shared_memory_bytes = 0;
  cudaStream_t const cuda_stream = nullptr;
  parallel::cuda_for_each<<<cuda_grid, cuda_block, shared_memory_bytes, cuda_stream>>>(f, first, last);
  
  // Explicit synchronization
  cudaDeviceSynchronize();
  
#elif defined(USE_OPENMP)
  // CPU execution path with OpenMP
  #pragma omp parallel for
  for (auto i = 0; i < n; ++i) {
    f(first[i]);
  }
#else
  // Sequential fallback
  for (; first != last; ++first) {
    f(*first);
  }
#endif
}

} // namespace parallel

struct AssignFunctor {
    double* output;
    double* input;

    // Constructor
    __host__ __device__ AssignFunctor(double* out, double* in)
        : output(out), input(in) {}

    // Overloaded function call operator
    __host__ __device__ void operator()(double& x) const {
        int index = &x - input;  // Compute index based on pointer arithmetic
        output[index] = x;       // Assign value
    }
};

struct AssignValueFunctor {
    double value;  // The value to assign

    // Constructor to initialize the value to assign
    __host__ __device__
    AssignValueFunctor(double v) : value(v) {}

    // Function call operator for the functor
    __host__ __device__
    void operator()(double& x) const {
        x = value;  // Assign the value to the array element
    }
};


#endif // PARALLEL_FOR_EACH_H

