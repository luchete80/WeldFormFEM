
#ifndef PARALLEL_FOR_EACH_H
#define PARALLEL_FOR_EACH_H

#include <iterator>
#include <utility>

// Preprocessor flags for enabling CUDA or OpenMP
#if defined(USE_CUDA)
#include <cuda_runtime.h>
#endif

#if defined(USE_OPENMP)
#include <omp.h>
#endif

namespace parallel {

#if defined(USE_CUDA)

/*
// CUDA kernel for parallel for_each
template <class F, class ForwardIt>
__global__
void cuda_for_each(F f, ForwardIt first, ForwardIt last) {
  using difference_type = typename std::iterator_traits<ForwardIt>::difference_type;
  auto const i = static_cast<difference_type>(threadIdx.x + blockIdx.x * blockDim.x);
  ForwardIt const it = first + i;
  if (it < last) f(*it);
}
*/

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
  if (n <= 0) return;

#if defined(USE_CUDA)
  // GPU execution path
  dim3 const cuda_block(32, 1, 1);
  dim3 const cuda_grid(parallel::ceildiv(static_cast<unsigned>(n), cuda_block.x), 1, 1);
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

#endif // PARALLEL_FOR_EACH_H




/////////////////////////////////////////////////////////
//////////////////////// FOR CONFIGURABLE BLOCK SIZE
/////////////////////////////////////////////////////////
/*
#ifndef PARALLEL_FOR_EACH_H
#define PARALLEL_FOR_EACH_H

#include <iterator>
#include <utility>

// Preprocessor flags for enabling CUDA or OpenMP
#if defined(USE_CUDA)
#include <cuda_runtime.h>
#endif

#if defined(USE_OPENMP)
#include <omp.h>
#endif

namespace parallel {

#if defined(USE_CUDA)

// CUDA kernel for parallel for_each
template <class F, class ForwardIt>
__global__
void cuda_for_each(F f, ForwardIt first, ForwardIt last);

// Helper function for ceiling division
template <class T>
__host__ __device__ inline constexpr
T ceildiv(T a, T b);

#endif

// Unified for_each function with configurable block size
template <class ForwardIt, class UnaryFunction>
void for_each(ForwardIt first, ForwardIt last, UnaryFunction f, unsigned block_size = 32);

} // namespace parallel

#endif // PARALLEL_FOR_EACH_H
*/


///////////////////////////////////////////////////////////////////
/*
//// POSSIOBLY FIXING FOR DEV FUNCTORS
#ifndef PARALLEL_FOR_EACH_H
#define PARALLEL_FOR_EACH_H

#include <iterator>
#include <utility>

// Preprocessor flags for enabling CUDA or OpenMP
#if defined(USE_CUDA)
#include <cuda_runtime.h>
#endif

#if defined(USE_OPENMP)
#include <omp.h>
#endif

namespace parallel {

#if defined(USE_CUDA)

// Explicitly marking the functor as `__device__`
template <typename F>
struct DeviceWrapper {
    F f;
    __device__ void operator()(double& x) { f(x); }
};

// CUDA kernel for `parallel_for_each` specialized for `double*`
template <class F>
__global__
void cuda_for_each(F f, double* first, double* last) {
    auto const i = threadIdx.x + blockIdx.x * blockDim.x;
    if (first + i < last) {
        f(first[i]);  // Ensuring `f` is a `__device__` function
    }
}

// Helper function for ceiling division
template <class T>
__host__ __device__ inline constexpr
T ceildiv(T a, T b) {
    return (a / b) + ((a % b) ? 1 : 0);
}

#endif // USE_CUDA

// Unified `for_each` function
template <class ForwardIt, class UnaryFunction>
void for_each(ForwardIt first, ForwardIt last, UnaryFunction f) {
    auto const n = last - first;
    if (n <= 0) return;

#if defined(USE_CUDA)
    // GPU execution path
    dim3 const cuda_block(32, 1, 1);
    dim3 const cuda_grid(parallel::ceildiv(static_cast<unsigned>(n), cuda_block.x), 1, 1);
    
    // Wrap the functor for device execution
    DeviceWrapper<UnaryFunction> d_f{f};

    cuda_for_each<<<cuda_grid, cuda_block>>>(d_f, first, last);
    
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

#endif // PARALLEL_FOR_EACH_H

*/
